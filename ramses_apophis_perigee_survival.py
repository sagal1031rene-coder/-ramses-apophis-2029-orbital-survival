# ==============================================================================
# SPACECRAFT RENDEZVOUS DYNAMICS: SOLAR RADIATION PRESSURE (SRP) & CRITICAL SWEEP
# Target: ESA Ramses Mission at Asteroid (99942) Apophis (April 2029 Perigee)
# Author: René Sagal Andrade (ORCID: 0009-0003-5574-5384)
# ==============================================================================

import rebound
import numpy as np
import matplotlib.pyplot as plt
from astroquery.jplhorizons import Horizons

def run_srp_critical_sweep():
    print("[INIT] Fetching April 2029 state vectors from NASA JPL Horizons...")
    epoch_jd = 2462236.5  # April 10, 2029
    
    # Cache solar system bodies to accelerate the critical altitude sweep
    bodies_info = []
    masses = {'10': 1.0, '399': 3.003489e-6, '301': 3.694e-8}
    for b_id, m in masses.items():
        obj = Horizons(id=b_id, location='@0', epochs=epoch_jd).vectors()
        bodies_info.append((m, float(obj['x']), float(obj['y']), float(obj['z']),
                            float(obj['vx']), float(obj['vy']), float(obj['vz'])))
        
    apophis = Horizons(id='99942', location='@0', epochs=epoch_jd).vectors()
    ax0, ay0, az0 = float(apophis['x']), float(apophis['y']), float(apophis['z'])
    avx0, avy0, avz0 = float(apophis['vx']), float(apophis['vy']), float(apophis['vz'])
    
    M_APOPHIS = 6.1e10 / 1.98847e30  # ~3.067e-20 Msun
    AU_TO_KM = 149597870.7
    AU_TO_M = 1.495978707e11
    DAY_TO_SEC = 86400.0
    
    # Spacecraft physical parameters (ESA Ramses class)
    # P_1AU = 4.56e-6 N/m^2 (Solar radiation pressure at 1 AU)
    # C_R = 1.3 (Reflectivity coefficient), A/m = 0.01 m^2/kg (10 m^2 / 1000 kg)
    P_1AU = 4.56e-6
    C_R = 1.3
    AREA_TO_MASS = 0.01  # m^2/kg
    
    # Convert SRP acceleration at 1 AU from m/s^2 to AU/day^2
    a_srp_1au_ms2 = P_1AU * C_R * AREA_TO_MASS
    a_srp_1au_auday2 = a_srp_1au_ms2 * (DAY_TO_SEC**2) / AU_TO_M
    
    # Sweep initial station-keeping radii from 0.6 km to 4.4 km (step 0.2 km)
    radii_sweep_km = np.arange(0.6, 4.6, 0.2)
    max_excursion_no_srp = []
    max_excursion_with_srp = []
    
    times = np.linspace(0, 7.0, 2500)
    
    print(f"[SWEEP] Running {len(radii_sweep_km)} altitude scenarios (Gravity-Only vs. SRP)...")
    for r0_km in radii_sweep_km:
        for use_srp in [False, True]:
            sim = rebound.Simulation()
            sim.units = ('AU', 'days', 'Msun')
            sim.integrator = "ias15"
            sim.ri_ias15.epsilon = 1e-12
            
            # Add Sun (0), Earth (1), Moon (2)
            for m, x, y, z, vx, vy, vz in bodies_info:
                sim.add(m=m, x=x, y=y, z=z, vx=vx, vy=vy, vz=vz)
                
            # Add Apophis (3)
            sim.add(m=M_APOPHIS, x=ax0, y=ay0, z=az0, vx=avx0, vy=avy0, vz=avz0)
            
            # Add Ramses spacecraft (4) in circular terminator-plane orbit (XZ plane)
            r0_au = r0_km / AU_TO_KM
            v_circ = np.sqrt(sim.G * M_APOPHIS / r0_au)
            sim.add(m=0.0, x=ax0 + r0_au, y=ay0, z=az0,
                    vx=avx0, vy=avy0, vz=avz0 + v_circ)
            
            # Inject Solar Radiation Pressure force on spacecraft (index 4)
            if use_srp:
                def srp_force(reb_sim):
                    particles = reb_sim.particles
                    sun = particles[0]
                    sc = particles[4]
                    dx = sc.x - sun.x
                    dy = sc.y - sun.y
                    dz = sc.z - sun.z
                    r_sun_au = np.sqrt(dx**2 + dy**2 + dz**2)
                    # Inverse-square solar flux scaling
                    acc_mag = a_srp_1au_auday2 / (r_sun_au**3)
                    sc.ax += acc_mag * dx
                    sc.ay += acc_mag * dy
                    sc.az += acc_mag * dz
                sim.additional_forces = srp_force
            
            max_dist_km = 0.0
            for t in times:
                sim.integrate(t)
                ast = sim.particles
                sc = sim.particles
                dist_km = np.sqrt((sc.x - ast.x)**2 + (sc.y - ast.y)**2 + (sc.z - ast.z)**2) * AU_TO_KM
                if dist_km > max_dist_km:
                    max_dist_km = dist_km
                    
            if use_srp:
                max_excursion_with_srp.append(max_dist_km)
            else:
                max_excursion_no_srp.append(max_dist_km)

    max_excursion_no_srp = np.array(max_excursion_no_srp)
    max_excursion_with_srp = np.array(max_excursion_with_srp)
    
    # Define survival threshold: spacecraft stays within 6.0 km bound regime post-flyby
    bound_no_srp = radii_sweep_km[max_excursion_no_srp < 6.0]
    bound_with_srp = radii_sweep_km[max_excursion_with_srp < 6.0]
    
    r_crit_no_srp = np.max(bound_no_srp) if len(bound_no_srp) > 0 else 0.0
    r_crit_with_srp = np.max(bound_with_srp) if len(bound_with_srp) > 0 else 0.0
    
    print("\n--- CRITICAL SURVIVAL ALTITUDE DIAGNOSTICS ---")
    print(f"SRP Acceleration at 1 AU      : {a_srp_1au_ms2:.3e} m/s^2")
    print(f"Max Safe Radius (Gravity-Only): {r_crit_no_srp:.2f} km")
    print(f"Max Safe Radius (With SRP)    : {r_crit_with_srp:.2f} km")
    print(f"SRP Safety Margin Reduction   : {r_crit_no_srp - r_crit_with_srp:.2f} km")

    # Generate Figure 2: Critical Altitude Sweep & SRP Perturbation Profile
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.5, 5.4), dpi=300)
    
    # Panel A: Maximum Post-Flyby Excursion vs. Initial Orbit Radius
    ax1.plot(radii_sweep_km, max_excursion_no_srp, color='#2980B9', marker='o', markersize=4.5,
             linewidth=2, label='Gravity-Only Baseline')
    ax1.plot(radii_sweep_km, max_excursion_with_srp, color='#D35400', marker='s', markersize=4.5,
             linewidth=2, linestyle='--', label='Coupled Gravity + Solar Radiation Pressure (SRP)')
    ax1.axhline(6.0, color='#E74C3C', linestyle=':', linewidth=1.5, label='Tidal Stripping Escape Threshold (6.0 km)')
    
    ax1.set_yscale('log')
    ax1.set_ylim(0.5, 500)
    ax1.set_title('A) Critical Survival Radius: Gravity vs. SRP Perturbation', fontweight='bold', fontsize=10.5)
    ax1.set_xlabel('Initial Station-Keeping Orbit Radius $r_0$ [km]')
    ax1.set_ylabel('Maximum Apophis-Centric Excursion [km] (Log)')
    ax1.grid(True, which='both', linestyle=':', alpha=0.55)
    ax1.legend(fontsize=8, loc='upper left')
    
    # Panel B: Acceleration Ratio (Solar Radiation Pressure vs. Apophis Gravity)
    r_profile = np.linspace(0.3, 5.0, 200)
    G_SI = 6.67430e-11
    M_AST_KG = 6.1e10
    a_grav_ast = (G_SI * M_AST_KG) / ((r_profile * 1000.0)**2)
    
    ax2.plot(r_profile, a_grav_ast * 1e6, color='#27AE60', linewidth=2.2, label=r"Apophis Local Gravity ($g_{\mathrm{ast}}$)")
    ax2.axhline(a_srp_1au_ms2 * 1e6, color='#F39C12', linestyle='--', linewidth=2.0,
                label=rf"Ramses Solar Panel SRP ({a_srp_1au_ms2*1e6:.3f} $\mu\mathrm{{m/s}}^2$)")
    
    # Mark equilibrium altitude where a_grav == a_srp
    r_eq_km = np.sqrt((G_SI * M_AST_KG) / a_srp_1au_ms2) / 1000.0
    ax2.axvline(r_eq_km, color='#8E44AD', linestyle=':', linewidth=1.6,    particles = sim.particles
    ast = particles[3]
    sc = particles[4]
            label=f'SRP Dominance Boundary ($r_{{eq}} = {r_eq_km:.2f}$ km)')
    
    ax2.set_yscale('log')
    ax2.set_title('B) Non-Gravitational SRP vs. Asteroid Gravity Balance', fontweight='bold', fontsize=10.5)
    ax2.set_xlabel('Radial Distance from Apophis Center [km]')
    ax2.set_ylabel(r'Acceleration Magnitude [$\mu\mathrm{m/s}^2$] (Log)')
    ax2.grid(True, which='both', linestyle=':', alpha=0.55)
    ax2.legend(fontsize=8, loc='upper right')
    
    plt.figtext(0.5, -0.02,
                'Source: REBOUND/IAS15 SRP & Tidal Sweep — René Sagal Andrade (ORCID: 0009-0003-5574-5384)',
                ha='center', fontsize=8, style='italic', color='#566573')
    
    plt.tight_layout()
    plt.savefig('ramses_srp_critical_sweep.png', bbox_inches='tight', dpi=300)
    print("\n[SUCCESS] Saved critical sweep chart as 'ramses_srp_critical_sweep.png'")
    plt.show()

if __name__ == "__main__":
    run_srp_critical_sweep()
