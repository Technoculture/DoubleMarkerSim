from down_syndrome_screening.models import create_default_population
import numpy as np


def summary_statistics(data: np.ndarray, verbose: bool = True):
    mean = np.mean(data)
    median = np.median(data)
    min_val = np.min(data)
    max_val = np.max(data)

    if verbose:
        print(f"\tMean: {mean:.2f}")
        print(f"\tMedian: {median:.2f}")
        print(f"\tMin: {min_val:.2f}")
        print(f"\tMax: {max_val:.2f}")

    return mean, median, min_val, max_val


# Create the default population
population = create_default_population()

# Generate a sample of 10,000 individuals
sample_size = 10_000
maternal_ages, marker_values, is_down = population.generate_sample(sample_size)

print(f"{np.sum(is_down)} Down syndrome cases out of {sample_size} total")
print(f"{maternal_ages[:5]=}")
print(f"{marker_values[:5]=}")
print(f"{is_down[:5]=}")

# Example: Calculate the mean age of mothers with Down syndrome pregnancies
print("Maternal age for Down syndrome pregnancies:")
down = summary_statistics(maternal_ages[is_down])
print("Maternal age for control pregnancies:")
control = summary_statistics(maternal_ages[~is_down])

print(f"Mean maternal age difference: {down[0] - control[0]:.2f}")
print(f"Median maternal age difference: {down[1] - control[1]:.2f}")

# Example: Calculate the median PAPP-A MoM for Down syndrome and control pregnancies
papp_a_index = 0  # Assuming PAPP-A is the first marker
print("PAPP-A MoM for Down syndrome pregnancies:")
pappa_down_stats = summary_statistics(np.exp(marker_values[is_down, papp_a_index]))

print("PAPP-A MoM for control pregnancies:")
papp_a_control = summary_statistics(np.exp(marker_values[~is_down, papp_a_index]))

# Example: Calculate correlation between PAPP-A and Free β-hCG for Down syndrome cases
corr_papp_a_bhcg_down = np.corrcoef(
    marker_values[is_down, 0], marker_values[is_down, 1]
)[0, 1]
print(f"Correlation between PAPP-A and Free β-hCG in Down syndrome cases: {
      corr_papp_a_bhcg_down:.3f}")
