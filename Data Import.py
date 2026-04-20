import pandas as pd
from astroquery.utils.tap.core import TapPlus
from astropy.coordinates import SkyCoord
import astropy.units as u
import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u

import csv

from scipy.__config__ import show
inputval = input('Show plots? (Y/N): ').lower()
if inputval == 'y':
    show_plots = True
else:
    show_plots = False
inputval = input('Path for saving data (None is default):')
if inputval.lower() == 'none' or inputval.lower() == 'data':
    save_path = "Data\\"
else:
    save_path = inputval
inputval = input('Save plots? (Y/N): ').lower()
if inputval == 'y':
    save_plots = True
else:
    save_plots = False
inputval = input('File containing cluster data:')
if inputval.lower() == 'none':
    cluster_file = "clusters_gaia_hms.csv"
else:
    cluster_file = inputval
clusters = []
afstanden = []
ra_list = []
dec_list = []
stralen = []

with open(cluster_file, newline="") as f:
    reader = csv.DictReader(f)
    for row in reader:
        clusters.append(row["Cluster"])
        afstanden.append(float(row["Afstand_pc"]))
        ra_list.append(row["RA_hms"])
        dec_list.append(row["DEC_dms"])
        stralen.append(float(row["Straal_arcmin"]))

print(clusters)
print(afstanden)
print(ra_list)
print(dec_list)
print(stralen)

for index, cluster in enumerate(clusters):

    print(f"Importing cluster: {cluster}")
    print(f"Right Ascension: {ra_list[index]}")
    # %% Parameters
    # coord = SkyCoord(ra=ra_list[index]*u.degree, dec=dec_list[index]*u.degree, frame="icrs")
    coord = SkyCoord(ra=ra_list[index], dec=dec_list[index], frame="icrs")
    ra0, dec0 = coord.ra.deg, coord.dec.deg
    rdeg      = stralen[index] / 60       # 20 arcmin naar graden
    d0_pc     = afstanden[index]          # clusterafstand (pc)
    radius_pc = d0_pc * np.radians(rdeg) 
    # radius_pc = 16

    # %% ADQL query
    ADQL = f"""
    SELECT source_id, ra, dec, parallax, parallax_error,
        pmra, pmra_error, pmdec, pmdec_error,
        radial_velocity, radial_velocity_error,
        rv_nb_transits, rv_chisq_pvalue, rv_renormalised_gof,
        phot_g_mean_mag, bp_rp,
        ruwe
    FROM gaiadr3.gaia_source
    WHERE 1 = CONTAINS(POINT('ICRS', ra, dec), CIRCLE('ICRS', {ra0}, {dec0}, {rdeg}))
    AND parallax > 0
    """

    tap = TapPlus(url="https://gea.esac.esa.int/tap-server/tap")
    df  = tap.launch_job_async(ADQL, verbose=False).get_results().to_pandas()
    df.columns = [c.lower() for c in df.columns]
    df["distance_pc"] = 1000.0 / df["parallax"]

    # %% Filteroverzicht (vóór 3D-bol en RUWE-filter)
    filters = {
        "Total in cone":                len(df),
        # "ruwe < 1.6":                   (df["ruwe"] < 20).sum(),
        # "radial_velocity not null":     df["radial_velocity"].notna().sum(),
        # "rv_nb_transits >= 10":         (df["rv_nb_transits"] >= 10).sum(),
        # "rv_chisq_pvalue not null":     df["rv_chisq_pvalue"].notna().sum(),
        # "rv_renormalised_gof not null": df["rv_renormalised_gof"].notna().sum(),
    }
    for label, count in filters.items():
        print(f"{label:<35} {count}")

    # %% RUWE filter
    # RUWE_MAX = 20
    # df_clean = df[df["ruwe"] < RUWE_MAX].copy()
    # print(f"Na RUWE < {RUWE_MAX} filter:         {len(df_clean)}  "
    #       f"(verwijderd: {len(df) - len(df_clean)})")

    # %% 3D bolfilter
    center      = SkyCoord(ra=ra0*u.deg, dec=dec0*u.deg, distance=d0_pc*u.pc, frame="icrs")
    stars       = SkyCoord(ra=df["ra"].values*u.deg, dec=df["dec"].values*u.deg,
                        distance=df["distance_pc"].values*u.pc, frame="icrs")
    distance_3d = stars.separation_3d(center).pc
    df = df[distance_3d < radius_pc].copy()
    print(f"\nStars in 3D sphere:          {len(df)}")

    df_clean = df.copy()

    # %% Sla beide datasets op
    df.to_csv(save_path + "data " + cluster + ".csv", index=False)
    print(f"Saved: data {cluster}.csv ({len(df_clean)} sterren)")

    # %% Sla variabelen op (gebruik df_clean als werkdataset)
    ra          = df_clean['ra']
    dec         = df_clean['dec']
    parallax    = df_clean['parallax']
    parallax_err = df_clean['parallax_error']
    pmra        = df_clean['pmra']
    pmra_err    = df_clean['pmra_error']
    pmdec       = df_clean['pmdec']
    pmdec_err   = df_clean['pmdec_error']
    rv          = df_clean['radial_velocity']
    rv_err      = df_clean['radial_velocity_error']
    rv_nb_transits  = df_clean['rv_nb_transits']
    rv_chisq_pvalue = df_clean['rv_chisq_pvalue']
    rv_renorm_gof   = df_clean['rv_renormalised_gof']
    g_mag       = df_clean['phot_g_mean_mag']
    bp_rp       = df_clean['bp_rp']
    ruwe        = df_clean['ruwe']

    # %% Sky map gekleurd op radiale snelheid
    if show_plots:
        fig, ax = plt.subplots(figsize=(7, 6))
        sc = ax.scatter(df_clean["ra"], df_clean["dec"],
                        c=df_clean["radial_velocity"], cmap="coolwarm",
                        s=8, alpha=0.8, linewidths=0)
        ax.invert_xaxis()
        ax.set_xlabel("RA [deg]")
        ax.set_ylabel("Dec [deg]")
        ax.set_title(f"Map of imported stars of {cluster}\n(RUWE < 1.6)")
        fig.colorbar(sc, ax=ax, label="Radial velocity [km s$^{-1}$]")
        plt.tight_layout()
        if save_plots:
            plt.savefig(save_path + f"Figures\\SkyPlots\\{cluster}_skyplot_ruwe.png", dpi=150)
            print(f"Saved: {cluster}_skyplot_ruwe.png")
        plt.show()

    # %% RUWE verdeling
    if show_plots:
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.hist(df["ruwe"].dropna(), bins=60, color="steelblue", edgecolor="none", alpha=0.8)
        # ax.axvline(RUWE_MAX, color="red", lw=1.5, ls="--", label=f"RUWE = {RUWE_MAX}")
        ax.set_xlabel("RUWE")
        ax.set_ylabel("Number of Stars")
        ax.set_title(f"RUWE-distribution — {cluster} (3D sphere)")
        ax.legend()
        plt.tight_layout()
        if save_plots:
            plt.savefig(save_path + f"Figures\\RUWEImportPlots\\{cluster}_ruwe_hist.png", dpi=150)
            print(f"Saved: {cluster}_ruwe_hist.png")
        plt.show()

print("\nFinished importing all clusters.")