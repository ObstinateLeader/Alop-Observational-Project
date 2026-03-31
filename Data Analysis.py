import numpy as np
import matplotlib.pyplot as plt


# from ALOP_Observation_Project_Data_Import_Script import ra, dec, parallax, parallax_err, pmra, pmra_err, pmdec, pmdec_err, rv, rv_err

# oude import
# import the data from the csv file saved by the import script
# data = np.genfromtxt("ngc6366_gaia_dr3.csv", delimiter=",", names=True)
# ra, dec, parallax, parallax_err, pmra, pmra_err, pmdec, pmdec_err, rv, rv_err = (
#     data["ra"],
#     data["dec"],
#     data["parallax"],
#     data["parallax_error"],
#     data["pmra"],
#     data["pmra_error"],
#     data["pmdec"],
#     data["pmdec_error"],
#     data["radial_velocity"],
#     data["radial_velocity_error"],
# )
inputval = input('Show plots? (Y/N): ').lower()
if inputval == 'y':
    show_plots = True
else:
    show_plots = False

inputval = input('Path for saved data (None is default):')
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
# nieuwe import
data = np.genfromtxt(save_path + "ngc6366_gaia_dr3_ruwe.csv", delimiter=",", names=True)


def analysis(data, cluster_name):
    global save_plots, show_plots
    print(f"Executing analysis of cluster: {cluster_name}")
    ra, dec, parallax, parallax_err, pmra, pmra_err, pmdec, pmdec_err, rv, rv_err, rv_nb_transits, rv_chisq_pvalue, rv_renorm_gof, g_mag, bp_rp, ruwe = (
        data["ra"],
        data["dec"],
        data["parallax"],
        data["parallax_error"],
        data["pmra"],
        data["pmra_error"],
        data["pmdec"],
        data["pmdec_error"],
        data["radial_velocity"],
        data["radial_velocity_error"],
        data["rv_nb_transits"],
        data["rv_chisq_pvalue"],
        data["rv_renormalised_gof"],
        data["phot_g_mean_mag"],
        data["bp_rp"],
        data["ruwe"],

    )

    # Checking the raw data
    print(f"Number of stars: {len(ra)}")
    print(f"Radial Velocities: {rv}")
    print(f"RUWE values: {ruwe}")

    rv_original = rv
    # masking rv, de nan's eruit filteren
    valid_rv_mask = ~np.isnan(rv)
    rv = rv[valid_rv_mask]
    rv_err = rv_err[valid_rv_mask]

    # Calculating the mean radial velocity and its standard deviation
    mean_rv = np.mean(rv)
    std_rv = np.std(rv)

    # Plotting

    star_indices = range(len(rv))
    if show_plots:
        plt.errorbar(star_indices, rv, yerr=rv_err, fmt='o', capsize=5, label='Radial Velocity with Errors')
        plt.axhline(mean_rv, color='r', linestyle='--', label=f'Mean RV = {mean_rv:.2f} km/s')
        # plt.fill_between(star_indices, mean_rv - std_rv, mean_rv + std_rv, color='r', alpha=0.2, label=f'±1σ = {std_rv:.2f} km/s')
        plt.fill_between(star_indices, mean_rv - 2 * std_rv, mean_rv + 2 * std_rv, color='r', alpha=0.1,
                        label=f'±2σ = {2 * std_rv:.2f} km/s')
        plt.xlabel('Star Index')
        plt.ylabel('Radial Velocity [km/s]')
        plt.title(f'Radial Velocities of Stars in {cluster_name}')
        plt.legend()
        if save_plots:
            plt.savefig(save_path + f"Figures\\RVPlots\\{cluster_name}_rv_plot.png", dpi=150)
            print(f"Saved: {cluster_name}_rv_plot.png")
        plt.show()

        plt.scatter(star_indices, rv_err, color='orange', label='Radial Velocity Errors')
        plt.axhline(np.mean(rv_err), color='g', linestyle='--', label=f'Mean RV Error = {np.mean(rv_err):.2f} km/s')
        # plt.axhline(std_rv, color='b', linestyle='--', label=f'Standard Deviation of RV = {std_rv:.2f} km/s')
        plt.xlabel('Star Index')
        plt.ylabel('Radial Velocity Error [km/s]')
        plt.title(f'Radial Velocity Errors of Stars in {cluster_name}')
        plt.legend()
        if save_plots:
            plt.savefig(save_path + f"Figures\\RVPlots\\{cluster_name}_rv_error_plot.png", dpi=150)
            print(f"Saved: {cluster_name}_rv_error_plot.png")
        plt.show()

    # Star_indices with RV flagged binaries:
    valid_indices = np.where(~np.isnan(rv_original))[0]

    # 2. Work with the valid data
    rv_clean = rv_original[valid_indices]
    mean_rv = np.mean(rv_clean)
    std_rv = np.std(rv_clean)

    # 3. Create a mask for the binaries within this 'clean' subset
    binary_mask = np.abs(rv_clean - mean_rv) > (2 * std_rv)  # Usually 3-sigma

    # 4. Map the 'subset index' back to the 'original index'
    # This picks the numbers from valid_indices where the binary_mask is True
    indices_RV_flagged = valid_indices[binary_mask]

    print(f"Original row numbers of RV-flagged binaries: {indices_RV_flagged}")

    # plt.scatter(star_indices, parallax, color='purple', label='Parallax')
    # plt.xlabel('Star Index')
    # plt.ylabel('Parallax [mas]')
    # plt.title(f'Parallaxes of Stars in {cluster_name} 6366')
    # plt.legend()
    # plt.show()

    print("---------Radial Velocity Analysis Results---------")
    print(f"Number of stars with radial velocity measurements: {len(rv)}")
    print(f"Mean Radial Velocity: {mean_rv:.2f} km/s")
    print(f"Standard Deviation of Radial Velocity: {std_rv:.2f} km/s")
    print(f"Number of flagged binaries with 1σ deviation: {np.sum(np.abs(rv - mean_rv) > std_rv)}")
    print(f"Number of flagged binaries with 2σ deviation: {np.sum(np.abs(rv - mean_rv) > 2 * std_rv)}")
    print("----------------------------------------------")

    ### RUWE
    RUWE_threshold = 1.22  # Aanpassen
    RUWE_max = 100

    mask = (ruwe > RUWE_threshold) & (ruwe < RUWE_max)

    # Filtering
    binary_candidates_RUWE = ruwe[mask]
    number_of_binaries_RUWE = len(binary_candidates_RUWE)

    if show_plots:
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.hist(ruwe, bins=60, color="steelblue", edgecolor="none", alpha=0.8)
        ax.axvline(RUWE_threshold, color="red", lw=1.5, ls="--", label=f"RUWE-threshold = {RUWE_threshold}")
        ax.set_xlabel("RUWE")
        ax.set_ylabel("Aantal sterren")
        ax.set_title(f"RUWE-verdeling — {cluster_name} (3D bol)")
        ax.set_xscale("log")
        ax.legend()
        plt.tight_layout()
        if save_plots:
            plt.savefig(save_path + f"Figures\\RUWEPlots\\{cluster_name}_ruwe_hist.png", dpi=150)
            print(f"Saved: {cluster_name}_ruwe_hist.png")
        plt.show()

    indices_RUWE_flagged = np.where(mask)[0]
    common_stars = np.intersect1d(indices_RV_flagged, indices_RUWE_flagged)

    print("---------RUWE Analysis Results---------")
    print(f"Number of binaries: {number_of_binaries_RUWE}")

    # indices flagged binaries RUWE
    # indices_RUWE_flagged = ruwe[mask]

    # Indices of flagged binaries from original list

    print(f"RUWE indices: {indices_RUWE_flagged}")
    print(f"RV indices: {indices_RV_flagged}")

    # Stars found by both methods
    print(f"Indices of stars that were flagged as binary by both methods: {common_stars}")
    print(f"Number of binaries flagged by both methods: {len(common_stars)}")
    print("----------------------------------------------")

    # Proper motions??
    return number_of_binaries_RUWE


# data = np.genfromtxt("ngc6366_gaia_dr3_ruwe_andere.csv", delimiter=",", names=True)
# cluster_name = "M45"
# analysis(data, cluster_name)


import csv

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
        # ra_list.append(float(row["RA_deg"]))
        # dec_list.append(float(row["DEC_deg"]))
        stralen.append(float(row["Straal_arcmin"]))

print(clusters)
print(afstanden)
# print(ra_list)
# print(dec_list)
print(stralen)

# importing data
# cluster_names = ["NGC6366, M45"]
cluster_names = clusters
number_of_binaries_ruwe_list = []
for cluster_name in cluster_names:
    data = np.genfromtxt(save_path + "data" + cluster_name + ".csv", delimiter=",", names=True)
    number_of_binaries_ruwe = analysis(data, cluster_name)
    number_of_binaries_ruwe_list.append(number_of_binaries_ruwe)

print("----------------Final Results------------------")
for index, cluster_name in enumerate(cluster_names):
    print(f"Number of flagged binaries ruwe method: {number_of_binaries_ruwe_list[index]} from cluster {cluster_name}")
print("-----------------------------------------------")


