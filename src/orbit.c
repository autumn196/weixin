#include <orbit.h>
#include <math.h>
#include <stdio.h>

/* Local constants using km units for consistency with the codebase */
#define MU_ORBIT 398600.4418  /* km³/s² */
#define EPSILON_ORBIT 1e-12

int orbit_elements_to_state(OrbitalElements *elements, StateVector *state) {
    if (!elements || !state) return -1;
    
    /* Convert orbital elements (in km and degrees) to position/velocity (in meters and m/s) */
    /* Input: a (km), e, i (deg), omega_big (deg), omega_small (deg), m0 (deg) */
    /* Output: position (m), velocity (m/s) */
    
    double a = elements->a * 1000.0;        /* Convert km to m */
    double e = elements->e;
    double i = elements->i * M_PI / 180.0;  /* Convert deg to rad */
    double Omega = elements->omega_big * M_PI / 180.0;
    double omega = elements->omega_small * M_PI / 180.0;
    double M = elements->m0 * M_PI / 180.0;
    
    double mu = MU_ORBIT * 1e9;  /* Convert km³/s² to m³/s² */
    
    /* Solve Kepler's equation: M = E - e*sin(E) */
    double E = M;  /* Initial guess */
    for (int iter = 0; iter < 50; iter++) {
        double f = E - e * sin(E) - M;
        double df = 1.0 - e * cos(E);
        if (fabs(df) < EPSILON_ORBIT) break;
        double dE = f / df;
        E = E - dE;
        if (fabs(dE) < EPSILON_ORBIT) break;
    }
    
    /* Compute true anomaly */
    double cos_E = cos(E);
    double sin_E = sin(E);
    double sqrt_1_minus_e2 = sqrt(1.0 - e * e);
    double nu = atan2(sqrt_1_minus_e2 * sin_E, cos_E - e);
    
    /* Compute distance */
    double r = a * (1.0 - e * cos_E);
    
    /* Position in orbital plane (perifocal coordinates) */
    double cos_nu = cos(nu);
    double sin_nu = sin(nu);
    double x_pqw = r * cos_nu;
    double y_pqw = r * sin_nu;
    
    /* Velocity in orbital plane */
    double p = a * (1.0 - e * e);
    double h = sqrt(mu * p);
    double vx_pqw = -(mu / h) * sin_nu;
    double vy_pqw = (mu / h) * (e + cos_nu);
    
    /* Rotation matrices components */
    double cos_Omega = cos(Omega);
    double sin_Omega = sin(Omega);
    double cos_omega = cos(omega);
    double sin_omega = sin(omega);
    double cos_i = cos(i);
    double sin_i = sin(i);
    
    /* Transform to ECI coordinates */
    /* R = Rz(-Omega) * Rx(-i) * Rz(-omega) */
    double r11 = cos_Omega * cos_omega - sin_Omega * sin_omega * cos_i;
    double r12 = -cos_Omega * sin_omega - sin_Omega * cos_omega * cos_i;
    double r21 = sin_Omega * cos_omega + cos_Omega * sin_omega * cos_i;
    double r22 = -sin_Omega * sin_omega + cos_Omega * cos_omega * cos_i;
    double r31 = sin_omega * sin_i;
    double r32 = cos_omega * sin_i;
    
    /* Apply rotation to position */
    state->position.x = r11 * x_pqw + r12 * y_pqw;
    state->position.y = r21 * x_pqw + r22 * y_pqw;
    state->position.z = r31 * x_pqw + r32 * y_pqw;
    
    /* Apply rotation to velocity */
    state->velocity.x = r11 * vx_pqw + r12 * vy_pqw;
    state->velocity.y = r21 * vx_pqw + r22 * vy_pqw;
    state->velocity.z = r31 * vx_pqw + r32 * vy_pqw;
    
    return 0;
}

int orbit_state_to_elements(StateVector *state, OrbitalElements *elements) {
    if (!state || !elements) return -1;
    
    /* Convert position/velocity (in meters and m/s) to orbital elements (in km and degrees) */
    double mu = MU_ORBIT * 1e9;  /* m³/s² */
    
    double x = state->position.x;
    double y = state->position.y;
    double z = state->position.z;
    double vx = state->velocity.x;
    double vy = state->velocity.y;
    double vz = state->velocity.z;
    
    double r_mag = sqrt(x*x + y*y + z*z);
    double v_mag = sqrt(vx*vx + vy*vy + vz*vz);
    
    if (r_mag < EPSILON_ORBIT) return -1;
    
    /* Angular momentum vector h = r × v */
    double hx = y * vz - z * vy;
    double hy = z * vx - x * vz;
    double hz = x * vy - y * vx;
    double h_mag = sqrt(hx*hx + hy*hy + hz*hz);
    
    if (h_mag < EPSILON_ORBIT) return -1;
    
    /* Node vector n = k × h (k = [0,0,1]) */
    double nx = -hy;
    double ny = hx;
    double n_mag = sqrt(nx*nx + ny*ny);
    
    /* Eccentricity vector e = (v × h)/μ - r/|r| */
    double r_dot_v = x*vx + y*vy + z*vz;
    double ex = (v_mag*v_mag - mu/r_mag) * x / mu - r_dot_v * vx / mu;
    double ey = (v_mag*v_mag - mu/r_mag) * y / mu - r_dot_v * vy / mu;
    double ez = (v_mag*v_mag - mu/r_mag) * z / mu - r_dot_v * vz / mu;
    double e = sqrt(ex*ex + ey*ey + ez*ez);
    
    /* Semi-major axis */
    double specific_energy = v_mag*v_mag/2.0 - mu/r_mag;
    double a;
    if (fabs(specific_energy) > EPSILON_ORBIT) {
        a = -mu / (2.0 * specific_energy);
    } else {
        a = r_mag;  /* Parabolic, approximate */
    }
    
    /* Inclination */
    double i = acos(hz / h_mag);
    
    /* Right ascension of ascending node (RAAN) */
    double Omega;
    if (n_mag > EPSILON_ORBIT) {
        Omega = acos(nx / n_mag);
        if (ny < 0) Omega = 2.0 * M_PI - Omega;
    } else {
        Omega = 0.0;  /* Equatorial orbit */
    }
    
    /* Argument of periapsis */
    double omega;
    if (n_mag > EPSILON_ORBIT && e > EPSILON_ORBIT) {
        double n_dot_e = nx*ex + ny*ey;
        omega = acos(n_dot_e / (n_mag * e));
        if (ez < 0) omega = 2.0 * M_PI - omega;
    } else {
        omega = 0.0;
    }
    
    /* True anomaly */
    double nu;
    if (e > EPSILON_ORBIT) {
        double e_dot_r = ex*x + ey*y + ez*z;
        nu = acos(e_dot_r / (e * r_mag));
        if (r_dot_v < 0) nu = 2.0 * M_PI - nu;
    } else {
        /* Circular orbit */
        if (n_mag > EPSILON_ORBIT) {
            double n_dot_r = nx*x + ny*y;
            nu = acos(n_dot_r / (n_mag * r_mag));
            if (z < 0) nu = 2.0 * M_PI - nu;
        } else {
            nu = atan2(y, x);
            if (nu < 0) nu += 2.0 * M_PI;
        }
    }
    
    /* Eccentric anomaly */
    double cos_nu = cos(nu);
    double E = atan2(sqrt(1.0 - e*e) * sin(nu), e + cos_nu);
    
    /* Mean anomaly */
    double M = E - e * sin(E);
    if (M < 0) M += 2.0 * M_PI;
    
    /* Store results (convert to km and degrees) */
    elements->a = a / 1000.0;  /* Convert m to km */
    elements->e = e;
    elements->i = i * 180.0 / M_PI;
    elements->omega_big = Omega * 180.0 / M_PI;
    elements->omega_small = omega * 180.0 / M_PI;
    elements->m0 = M * 180.0 / M_PI;
    
    return 0;
}

Vector3 orbit_eccentricity_vector(Vector3 r, Vector3 v) {
    return (Vector3){0, 0, 0};
}

double orbit_period(double a) {
    if (a <= 0) return 0;
    double mu = 3.986004418e14;
    return 2 * M_PI * sqrt(a * a * a / mu);
}

double orbit_orbital_energy(double a) {
    if (a <= 0) return 0;
    double mu = 3.986004418e14;
    return -mu / (2 * a);
}

double orbit_velocity_circular(double semi_major_axis) {
    if (semi_major_axis <= 0) return 0;
    double mu = 3.986004418e14;
    return sqrt(mu / semi_major_axis);
}

double orbit_kepler_equation_solve(double M, double e, double tolerance, int max_iterations) {
    double E = M;
    for (int i = 0; i < max_iterations; i++) {
        double f = E - e * sin(E) - M;
        double df = 1 - e * cos(E);
        if (fabs(df) < 1e-15) break;
        double E_new = E - f / df;
        if (fabs(E_new - E) < tolerance) return E_new;
        E = E_new;
    }
    return E;
}

double orbit_true_anomaly_from_mean(double M, double e) {
    double E = orbit_kepler_equation_solve(M, e, 1e-12, 20);
    return 2 * atan2(sqrt(1+e) * sin(E/2), sqrt(1-e) * cos(E/2));
}

int orbit_hohmann_transfer(double a1, double a2, HohmannTransfer *transfer) {
    if (!transfer || a1 <= 0 || a2 <= 0) return -1;
    
    double mu = 3.986004418e14;
    double v1 = sqrt(mu / a1);
    double v2 = sqrt(mu / a2);
    double a_transfer = (a1 + a2) / 2;
    double v_periapsis = sqrt(mu * (2/a1 - 1/a_transfer));
    double v_apoapsis = sqrt(mu * (2/a2 - 1/a_transfer));
    
    transfer->delta_v1 = fabs(v_periapsis - v1);
    transfer->delta_v2 = fabs(v2 - v_apoapsis);
    transfer->total_delta_v = transfer->delta_v1 + transfer->delta_v2;
    transfer->transfer_time = M_PI * sqrt(a_transfer * a_transfer * a_transfer / mu);
    
    return 0;
}

int orbit_lambert_solve_simple(Vector3 r1, Vector3 r2, double tof, double mu, Vector3 *v1, Vector3 *v2) {
    if (!v1 || !v2 || tof <= 0 || mu <= 0) return -1;
    
    /* Simple Lambert solver using Battin's method approximation */
    double r1_mag = sqrt(r1.x*r1.x + r1.y*r1.y + r1.z*r1.z);
    double r2_mag = sqrt(r2.x*r2.x + r2.y*r2.y + r2.z*r2.z);
    
    if (r1_mag < EPSILON_ORBIT || r2_mag < EPSILON_ORBIT) return -1;
    
    /* Compute the chord and related quantities */
    double c = sqrt((r2.x-r1.x)*(r2.x-r1.x) + (r2.y-r1.y)*(r2.y-r1.y) + (r2.z-r1.z)*(r2.z-r1.z));
    double s = (r1_mag + r2_mag + c) / 2.0;
    
    /* Determine transfer angle */
    double cos_dnu = (r1.x*r2.x + r1.y*r2.y + r1.z*r2.z) / (r1_mag * r2_mag);
    if (cos_dnu > 1.0) cos_dnu = 1.0;
    if (cos_dnu < -1.0) cos_dnu = -1.0;
    
    /* Cross product for direction */
    double cross_z = r1.x*r2.y - r1.y*r2.x;
    int direction = (cross_z >= 0) ? 1 : -1;
    double sin_dnu = sqrt(1.0 - cos_dnu*cos_dnu) * direction;
    
    /* Compute semi-major axis using simplified approach */
    /* For short transfers, approximate with a minimum energy ellipse */
    double a_min = s / 2.0;
    double tof_min = sqrt(2.0) / 3.0 * sqrt(a_min*a_min*a_min / mu) * (M_PI - (s - c) / a_min * sqrt((s - c) / (2.0 * a_min)));
    
    /* Adjust semi-major axis based on transfer time ratio */
    double ratio = tof / tof_min;
    double a;
    if (ratio > 1.0) {
        a = a_min * pow(ratio, 2.0/3.0);
    } else {
        a = a_min;
    }
    
    /* Compute parameter p */
    double p = a * (1.0 - (r1_mag - r2_mag) * (r1_mag - r2_mag) / (c * c));
    if (p < EPSILON_ORBIT) p = 0.5 * (r1_mag + r2_mag);
    
    /* Compute f and g functions */
    double f = 1.0 - r2_mag / p * (1.0 - cos_dnu);
    double g = r1_mag * r2_mag * sin_dnu / sqrt(mu * p);
    double g_dot = 1.0 - r1_mag / p * (1.0 - cos_dnu);
    
    /* Compute velocities */
    if (fabs(g) < EPSILON_ORBIT) return -1;
    
    v1->x = (r2.x - f * r1.x) / g;
    v1->y = (r2.y - f * r1.y) / g;
    v1->z = (r2.z - f * r1.z) / g;
    
    v2->x = g_dot * v1->x + (r2.x - f * r1.x) / g * (1.0 - g_dot);
    v2->y = g_dot * v1->y + (r2.y - f * r1.y) / g * (1.0 - g_dot);
    v2->z = g_dot * v1->z + (r2.z - f * r1.z) / g * (1.0 - g_dot);
    
    /* Recalculate v2 using the correct formula */
    double f_dot = sqrt(mu / p) * tan((acos(cos_dnu)) / 2.0) * ((1.0 - cos_dnu) / p - 1.0 / r1_mag - 1.0 / r2_mag);
    v2->x = f_dot * r1.x + g_dot * v1->x;
    v2->y = f_dot * r1.y + g_dot * v1->y;
    v2->z = f_dot * r1.z + g_dot * v1->z;
    
    return 0;
}

int orbit_lambert_solve(Vector3 r1, Vector3 r2, double tof, int num_revolutions, double mu, double max_iterations, double convergence, Vector3 *v1, Vector3 *v2) {
    if (!v1 || !v2) return -1;
    (void)num_revolutions;
    (void)max_iterations;
    (void)convergence;
    
    /* Use simple solver for now */
    return orbit_lambert_solve_simple(r1, r2, tof, mu, v1, v2);
}

/**
 * Compute gravitational acceleration at position r (in meters)
 * Returns acceleration in m/s²
 */
static Vector3 compute_gravity_acceleration(Vector3 r) {
    double mu = MU_ORBIT * 1e9;  /* m³/s² */
    double r_mag = sqrt(r.x*r.x + r.y*r.y + r.z*r.z);
    
    if (r_mag < EPSILON_ORBIT) {
        return (Vector3){0, 0, 0};
    }
    
    double factor = -mu / (r_mag * r_mag * r_mag);
    return (Vector3){factor * r.x, factor * r.y, factor * r.z};
}

int orbit_rk4_step(StateVector *state, double dt, Vector3 *external_acceleration) {
    if (!state) return -1;
    
    /* RK4 integration for orbital mechanics */
    /* State: [x, y, z, vx, vy, vz] */
    /* Derivative: [vx, vy, vz, ax, ay, az] */
    
    Vector3 r = state->position;
    Vector3 v = state->velocity;
    
    /* k1 */
    Vector3 a1 = compute_gravity_acceleration(r);
    if (external_acceleration) {
        a1.x += external_acceleration->x;
        a1.y += external_acceleration->y;
        a1.z += external_acceleration->z;
    }
    Vector3 k1_r = v;
    Vector3 k1_v = a1;
    
    /* k2 */
    Vector3 r2 = {r.x + 0.5*dt*k1_r.x, r.y + 0.5*dt*k1_r.y, r.z + 0.5*dt*k1_r.z};
    Vector3 v2 = {v.x + 0.5*dt*k1_v.x, v.y + 0.5*dt*k1_v.y, v.z + 0.5*dt*k1_v.z};
    Vector3 a2 = compute_gravity_acceleration(r2);
    if (external_acceleration) {
        a2.x += external_acceleration->x;
        a2.y += external_acceleration->y;
        a2.z += external_acceleration->z;
    }
    Vector3 k2_r = v2;
    Vector3 k2_v = a2;
    
    /* k3 */
    Vector3 r3 = {r.x + 0.5*dt*k2_r.x, r.y + 0.5*dt*k2_r.y, r.z + 0.5*dt*k2_r.z};
    Vector3 v3 = {v.x + 0.5*dt*k2_v.x, v.y + 0.5*dt*k2_v.y, v.z + 0.5*dt*k2_v.z};
    Vector3 a3 = compute_gravity_acceleration(r3);
    if (external_acceleration) {
        a3.x += external_acceleration->x;
        a3.y += external_acceleration->y;
        a3.z += external_acceleration->z;
    }
    Vector3 k3_r = v3;
    Vector3 k3_v = a3;
    
    /* k4 */
    Vector3 r4 = {r.x + dt*k3_r.x, r.y + dt*k3_r.y, r.z + dt*k3_r.z};
    Vector3 v4 = {v.x + dt*k3_v.x, v.y + dt*k3_v.y, v.z + dt*k3_v.z};
    Vector3 a4 = compute_gravity_acceleration(r4);
    if (external_acceleration) {
        a4.x += external_acceleration->x;
        a4.y += external_acceleration->y;
        a4.z += external_acceleration->z;
    }
    Vector3 k4_r = v4;
    Vector3 k4_v = a4;
    
    /* Update state */
    state->position.x = r.x + dt/6.0 * (k1_r.x + 2.0*k2_r.x + 2.0*k3_r.x + k4_r.x);
    state->position.y = r.y + dt/6.0 * (k1_r.y + 2.0*k2_r.y + 2.0*k3_r.y + k4_r.y);
    state->position.z = r.z + dt/6.0 * (k1_r.z + 2.0*k2_r.z + 2.0*k3_r.z + k4_r.z);
    
    state->velocity.x = v.x + dt/6.0 * (k1_v.x + 2.0*k2_v.x + 2.0*k3_v.x + k4_v.x);
    state->velocity.y = v.y + dt/6.0 * (k1_v.y + 2.0*k2_v.y + 2.0*k3_v.y + k4_v.y);
    state->velocity.z = v.z + dt/6.0 * (k1_v.z + 2.0*k2_v.z + 2.0*k3_v.z + k4_v.z);
    
    state->time += dt;
    
    return 0;
}

int orbit_propagate(StateVector *initial_state, double propagation_time, StateVector *final_state, double time_step) {
    if (!initial_state || !final_state || propagation_time < 0) return -1;
    
    /* Copy initial state to final state */
    *final_state = *initial_state;
    
    if (propagation_time < EPSILON_ORBIT) return 0;
    
    /* Use default time step if not specified */
    if (time_step <= 0) time_step = 1.0;  /* 1 second default */
    
    /* Propagate using RK4 */
    double elapsed = 0;
    while (elapsed < propagation_time) {
        double dt = time_step;
        if (elapsed + dt > propagation_time) {
            dt = propagation_time - elapsed;
        }
        
        if (orbit_rk4_step(final_state, dt, NULL) != 0) {
            return -1;
        }
        
        elapsed += dt;
    }
    
    return 0;
}

double orbit_delta_v_inclination_change(double a, double di) {
    double mu = 3.986004418e14;
    double v = sqrt(mu / a);
    return 2 * v * sin(di / 2);
}

double orbit_delta_v_raan_change(double a, double e, double i, double dOmega) {
    double mu = 3.986004418e14;
    double v = sqrt(mu / a);
    return 2 * v * sin(i) * sin(dOmega / 2);
}

double orbit_delta_v_eccentricity_change(double a, double de) {
    return 0;
}

double orbit_closest_approach_distance(OrbitalElements *orb1, OrbitalElements *orb2) {
    if (!orb1 || !orb2) return 1e10;
    return 1e10;
}

int orbit_orbits_intersect(OrbitalElements *orb1, OrbitalElements *orb2, double *intersection_distance) {
    if (!orb1 || !orb2 || !intersection_distance) return 0;
    *intersection_distance = 1e10;
    return 0;
}

double orbit_decay_rate(double altitude, double ballistic_coefficient) {
    return 0;
}

double orbit_maintenance_delta_v(double decay_rate, double maintenance_period) {
    return 0;
}

int orbit_calculate_relative_elements(OrbitalElements *chief_orb, OrbitalElements *deputy_orb, OrbitalElements *relative_orb) {
    if (!chief_orb || !deputy_orb || !relative_orb) return -1;
    return 0;
}

int orbit_absolute_from_relative(OrbitalElements *chief_orb, OrbitalElements *relative_orb, OrbitalElements *deputy_orb) {
    if (!chief_orb || !relative_orb || !deputy_orb) return -1;
    return 0;
}

int orbit_lambert_find_optimal_tof(Vector3 r1, Vector3 r2, double tof_min, double tof_max, double mu, double *optimal_tof) {
    if (!optimal_tof) return -1;
    *optimal_tof = (tof_min + tof_max) / 2;
    return 0;
}

double orbit_delta_v_magnitude(Vector3 v1, Vector3 v2) {
    Vector3 dv = vector3_sub(v2, v1);
    return vector3_magnitude(dv);
}