# src/down_syndrome_screening/analysis.py

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from scipy.stats import norm
from models import create_default_population


def calculate_lr(
    log_value, log_median_control, log_sd_control, log_median_down, log_sd_down
):
    """Calculate likelihood ratio for log-transformed values."""
    lr = norm.pdf(log_value, log_median_down, log_sd_down) / norm.pdf(
        log_value, log_median_control, log_sd_control
    )
    return np.clip(lr, 0, 1e6)  # Cap the max LR at 1e6


def maternal_age_risk(age):
    """Calculate risk of Down syndrome based on maternal age."""
    return 1 / (1000 * np.exp(-0.1 * (age - 20)))  # More realistic age-related risk


def generate_and_analyze_population(sample_size: int = 100_000):
    population = create_default_population()
    maternal_ages, marker_values, is_down = population.generate_sample(sample_size)

    df = pd.DataFrame(
        {
            "Maternal Age": maternal_ages,
            "Free β-hCG": marker_values[:, 0],
            "PAPP-A": marker_values[:, 1],
            "NT": marker_values[:, 2],
            "Is Down": is_down,
        }
    )

    # Calculate likelihood ratios for each marker
    for marker in ["Free β-hCG", "PAPP-A", "NT"]:
        log_values = np.log(df[marker])
        log_median_control = np.log(df[~df["Is Down"]][marker].median())
        log_sd_control = df[~df["Is Down"]][marker].apply(np.log).std()
        log_median_down = np.log(df[df["Is Down"]][marker].median())
        log_sd_down = df[df["Is Down"]][marker].apply(np.log).std()

        print(f"\n{marker}:")
        print(
            f"Control - Median: {np.exp(log_median_control)                                 :.4f}, Log SD: {log_sd_control:.4f}"
        )
        print(
            f"Down Syndrome - Median: {np.exp(log_median_down)                                       :.4f}, Log SD: {log_sd_down:.4f}"
        )

        df[f"LR {marker}"] = log_values.apply(
            calculate_lr,
            args=(log_median_control, log_sd_control, log_median_down, log_sd_down),
        )

    # Calculate combined likelihood ratio
    df["Combined LR"] = df["LR Free β-hCG"] * df["LR PAPP-A"] * df["LR NT"]

    # Calculate age-related risk
    df["Age Risk"] = df["Maternal Age"].apply(maternal_age_risk)

    # Calculate final risk
    df["Final Risk"] = df["Age Risk"] * df["Combined LR"] / 100  # Adjust the scaling

    # Determine screen positive based on a risk threshold (e.g., 1 in 250)
    risk_threshold = 1 / 250
    df["Screen Positive"] = df["Final Risk"] > risk_threshold

    # Calculate detection rate and false positive rate
    detection_rate = (df["Screen Positive"] & df["Is Down"]).sum() / df["Is Down"].sum()
    false_positive_rate = (df["Screen Positive"] & ~df["Is Down"]).sum() / (
        ~df["Is Down"]
    ).sum()

    print(f"\nDetection Rate: {detection_rate:.2%}")
    print(f"False Positive Rate: {false_positive_rate:.2%}")

    # Plot ROC curve
    thresholds = np.logspace(-6, -1, 100)
    tpr = []
    fpr = []
    for threshold in thresholds:
        df["Screen Positive"] = df["Final Risk"] > threshold
        tpr.append((df["Screen Positive"] & df["Is Down"]).sum() / df["Is Down"].sum())
        fpr.append(
            (df["Screen Positive"] & ~df["Is Down"]).sum() / (~df["Is Down"]).sum()
        )

    plt.figure(figsize=(10, 8))
    plt.plot(fpr, tpr)
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate (Detection Rate)")
    plt.title("ROC Curve for Down Syndrome Screening")
    plt.show()

    # Print some statistics about the likelihood ratios
    for column in ["LR Free β-hCG", "LR PAPP-A", "LR NT", "Combined LR"]:
        print(f"\n{column}:")
        print(f"Min: {df[column].min():.4f}")
        print(f"Max: {df[column].max():.4f}")
        print(f"Mean: {df[column].mean():.4f}")
        print(f"Median: {df[column].median():.4f}")


if __name__ == "__main__":
    generate_and_analyze_population()

