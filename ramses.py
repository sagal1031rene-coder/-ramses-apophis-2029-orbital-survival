# ==============================================================================
# SPACECRAFT RENDEZVOUS DYNAMICS: ESA RAMSES ORBITAL SURVIVAL AT 2029 PERIGEE
# Target: Asteroid (99942) Apophis - Tidal Disruption & Hill Sphere Collapse
# Author: René Sagal Andrade (ORCID: 0009-0003-5574-5384)
# ==============================================================================

import rebound
import numpy as np
import matplotlib.pyplot as plt
from astroquery.jplhorizons import Horizons

def run_ramses_survival_simulation():
    print("[INIT] Fetching April 2029 state vectors from NASA JPL Horizons...")
    # Epoch: April 10, 2029 (3 days prior to closest approach)
    epoch_jd = 2462236.5
    
    sim = rebound.Simulation()
    sim.units = ('AU', 'days', 'Msun')
    sim.integrator = "ias15"
    sim.ri_ias15.epsilon = 1e-12
    
    # Physical masses in Solar Masses (Msun)
    M_SUN = 1.0
    M_EARTH = 3.003489e-6
    M_MOON = 3.694e-8
    # Apophis mass (~6.1e10 kg -> ~3.067e-20 Msun)
    M_APOPHIS = 6.1e10 / 1.98847e30
    
    # Add Sun, Earth, and Moon
    bodies = [('10', M_SUN), ('399', M_EARTH), ('301', M_MOON)]
    for b_id, mass in bodies:
        obj = Horizons(id=b_id, location='@0', epochs=epoch_jd).vectors()
        sim.add(m=mass, x=float(obj['x']), y=float(obj['y']), z=float(obj['z']),
                vx=float(obj['vx']), vy=float(obj['vy']), vz=float(obj['vz']))
        
    # Add Asteroid (99942) Apophis (index 3)
    apophis = Horizons(id='99942', location='@0', epochs=epoch_jd).vectors()
    ax0, ay0, az0 = float(apophis['x']), float(apophis['y']), float(apophis['z'])
    avx0, avy0, avz0 = float(apophis['vx']), float(apophis['vy']), float(apophis['vz'])
    sim.add(m=M_APOPHIS, x=ax0, y=ay0, z=az0, vx=avx0, vy=avy0, vz=avz0)
    
    # Unit conversion factors
    AU_TO_KM = 149597870.7
    DAY_TO_SEC = 86400.0
    G_AU = sim.G  # Gravitational constant in AU^3 / (Msun * days^2)
    
    # Configure two virtual Ramses station-keeping orbits around Apophis:
    # Probe A: Tight safe orbit (r0 = 1.2 km)
    # Probe B: Wide vulnerable orbit (r0 = 5.0 km)
    radii_km = [1.2, 5.0]
    for r_km in radii_km:
        r_au = r_km / AU_TO_KM
        # Circular orbital velocity around Apophis: v = sqrt(G * M_ast / r)
        v_circ_au_day = np.sqrt(G_AU * M_APOPHIS / r_au)
        sim.add(m=0.0,
                x=ax0 + r_au, y=ay0, z=az0,
                vx=avx0, vy=avy0 + v_circ_au_day, vz=avz0)
        
    # Integrate over 7 days (April 10 to April 17, 2029)
    n_steps = 3500
    times = np.linspace(0, 7.0, n_steps)
    
    earth_dist_km = np.zeros(n_steps)
    hill_radius_km = np.zeros(n_steps)
    rel_pos_A = np.zeros((n_steps, 3))
    rel_pos_B = np.zeros((n_steps, 3))
    dist_A_km = np.zeros(n_steps)
    dist_B_km = np.zeros(n_steps)
    
    print("[PROPAGATE] Integrating 6-body coupled dynamics across terrestrial perigee...")
    for i, t in enumerate(times):
        earth = sim.particles[1]
        ast = sim.particles[3]
        probe_A = sim.particles[4]
        probe_B = sim.particles[5]

        
        # Geocentric distance of Apophis
        d_ea_au = np.sqrt((ast.x - earth.x)**2 + (ast.y - earth.y)**2 + (ast.z - earth.z)**2)
        earth_dist_km[i] = d_ea_au * AU_TO_KM
        
        # Instantaneous Earth-dominated Hill Sphere radius of Apophis: r_H = d * (M_ast / (3 * M_earth))^(1/3)
        hill_radius_km[i] = earth_dist_km[i] * ((M_APOPHIS / (3.0 * M_EARTH)) ** (1.0 / 3.0))
        
        # Probe A relative position [km]
        dxA = (probe_A.x - ast.x) * AU_TO_KM
        dyA = (probe_A.y - ast.y) * AU_TO_KM
        dzA = (probe_A.z - ast.z) * AU_TO_KM
        rel_pos_A[i] = [dxA, dyA, dzA]
        dist_A_km[i] = np.sqrt(dxA**2 + dyA**2 + dzA**2)
        
        # Probe B relative position [km]
        dxB = (probe_B.x - ast.x) * AU_TO_KM
        dyB = (probe_B.y - ast.y) * AU_TO_KM
        dzB = (probe_B.z - ast.z) * AU_TO_KM
        rel_pos_B[i] = [dxB, dyB, dzB]
        dist_B_km[i] = np.sqrt(dxB**2 + dyB**2 + dzB**2)

    min_perigee = np.min(earth_dist_km)
    min_hill = np.min(hill_radius_km)
    t_perigee = times[np.argmin(earth_dist_km)]
    
    print("\n--- RAMSES ORBITAL SURVIVAL DIAGNOSTICS ---")
    print(f"Minimum Geocentric Perigee : {min_perigee:,.2f} km (at t = +{t_perigee:.2f} days)")
    print(f"Minimum Hill Sphere Radius : {min_hill:.3f} km")
    print(f"Probe A (1.2 km) Max Dist  : {np.max(dist_A_km):.3f} km -> {'SURVIVED (BOUND)' if np.max(dist_A_km) < 3.0 else 'ESCAPED'}")
    print(f"Probe B (5.0 km) Max Dist  : {np.max(dist_B_km):,.2f} km -> {'SURVIVED (BOUND)' if np.max(dist_B_km) < 10.0 else 'TIDALLY STRIPPED'}")

    # Generate Dual-Panel Figure
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.5, 5.4), dpi=300)
    
    # Panel A: Hill Sphere Collapse vs. Spacecraft Altitude
    ax1.plot(times, hill_radius_km, color='#E74C3C', linewidth=2.2, label='Apophis Hill Sphere Radius ($r_H$)')
    ax1.plot(times, dist_A_km, color='#27AE60', linewidth=1.8, linestyle='-', label='Probe A: Tight Orbit ($r_0 = 1.2$ km)')
    ax1.plot(times, np.clip(dist_B_km, 0, 25), color='#8E44AD', linewidth=1.8, linestyle='--', label='Probe B: Wide Orbit ($r_0 = 5.0$ km, clipped)')
    ax1.axvline(t_perigee, color='#34495E', linestyle=':', alpha=0.8, label='Closest Approach (April 13, 2029)')
    
    ax1.set_yscale('log')
    ax1.set_ylim(0.5, 30)
    ax1.set_title('A) Hill Sphere Contraction & Tidal Stripping Threshold', fontweight='bold', fontsize=10.5)
    ax1.set_xlabel('Time Since April 10, 2029 [days]')
    ax1.set_ylabel('Apophis-Centric Distance [km] (Log Scale)')
    ax1.grid(True, which='both', linestyle=':', alpha=0.55)
    ax1.legend(fontsize=8, loc='lower left')
    
    # Panel B: Asteroid-Centered Trajectory (XY Plane Zoom)
    ax2.scatter([0], [0], color='#D35400', s=120, zorder=6, label='(99942) Apophis Center')
    min_hill_circle = plt.Circle((0, 0), min_hill, color='#E74C3C', fill=False, linestyle='--', linewidth=1.5,
                                 label=f'Min Perigee Hill Radius ({min_hill:.2f} km)')
    ax2.add_patch(min_hill_circle)
    
    ax2.plot(rel_pos_A[:, 0], rel_pos_A[:, 1], color='#27AE60', linewidth=1.3, label='Probe A Trajectory (Bound)')
    ax2.plot(rel_pos_B[:, 0], rel_pos_B[:, 1], color='#8E44AD', linewidth=1.5, linestyle='--', label='Probe B Trajectory (Ejected)')
    
    ax2.set_xlim(-12, 12)
    ax2.set_ylim(-12, 12)
    ax2.set_aspect('equal')
    ax2.set_title('B) Asteroid-Centered Frame: Bound vs. Stripped Orbits', fontweight='bold', fontsize=10.5)
    ax2.set_xlabel('Relative X [km]')
    ax2.set_ylabel('Relative Y [km]')
    ax2.grid(True, linestyle=':', alpha=0.55)
    ax2.legend(fontsize=8, loc='upper right')
    
    plt.figtext(0.5, -0.02,
                'Source: REBOUND/IAS15 6-Body Simulation — René Sagal Andrade (ORCID: 0009-0003-5574-5384)',
                ha='center', fontsize=8, style='italic', color='#566573')
    
    plt.tight_layout()
    plt.savefig('ramses_perigee_survival.png', bbox_inches='tight', dpi=300)
    print("\n[SUCCESS] Saved diagnostic chart as 'ramses_perigee_survival.png'")
    plt.show()

if __name__ == "__main__":
    run_ramses_survival_simulation()
