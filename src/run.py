from down_syndrome_screening.models import create_default_population
import numpy as np

# Create the default population
population = create_default_population()

# Generate a sample of 10,000 individuals
sample_size = 10000
maternal_ages, marker_values, is_down = population.generate_sample(sample_size)

# Example: Calculate the mean age of mothers with Down syndrome pregnancies
mean_age_down = np.mean(maternal_ages[is_down])
print(f"Mean maternal age for Down syndrome pregnancies: {mean_age_down:.2f}")

# Example: Calculate the median PAPP-A MoM for Down syndrome and control pregnancies
papp_a_index = 0  # Assuming PAPP-A is the first marker
median_papp_a_down = np.median(np.exp(marker_values[is_down, papp_a_index]))
median_papp_a_control = np.median(np.exp(marker_values[~is_down, papp_a_index]))
print(f"Median PAPP-A MoM for Down syndrome pregnancies: {median_papp_a_down:.2f}")
print(f"Median PAPP-A MoM for control pregnancies: {median_papp_a_control:.2f}")

# Example: Calculate correlation between PAPP-A and Free β-hCG for Down syndrome cases
corr_papp_a_bhcg_down = np.corrcoef(
    marker_values[is_down, 0], marker_values[is_down, 1]
)[0, 1]
print(f"Correlation between PAPP-A and Free β-hCG in Down syndrome cases: {
      corr_papp_a_bhcg_down:.3f}")

