import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from models import create_default_population


def generate_and_visualize_population(sample_size: int = 100_000):
    # Create the default population
    population = create_default_population()

    # Generate a sample
    maternal_ages, marker_values, is_down = population.generate_sample(sample_size)

    # Print diagnostic information
    print(f"Total samples: {sample_size}")
    print(f"Down syndrome cases: {np.sum(is_down)}")
    print(f"Control cases: {np.sum(~is_down)}")
    print(f"Observed prevalence: {np.mean(is_down):.6f}")

    # Create a DataFrame with log MoM values
    df = pd.DataFrame(
        {
            "Log Free β-hCG": marker_values[:, 1],
            "Log PAPP-A": marker_values[:, 0],
            "Log NT": marker_values[:, 6],
            "Group": ["Down Syndrome" if d else "Control" for d in is_down],
        }
    )

    # Set up the plot
    fig, axs = plt.subplots(3, 1, figsize=(12, 18))
    fig.suptitle(f"Distribution of Log MoM Markers in a Population of {
                 sample_size:,}")

    # Color palette
    colors = {"Control": "blue", "Down Syndrome": "red"}

    # Function to plot histogram
    def plot_histogram(ax, marker):
        # Determine the x-axis range for both groups combined
        all_data = df[marker]
        lower, upper = all_data.min(), all_data.max()
        print(f"{marker}: {lower=}, {upper=}")

        for group in ["Control", "Down Syndrome"]:
            data = df[df["Group"] == group][marker]
            sns.histplot(
                data=data,
                kde=True,
                ax=ax,
                color=colors[group],
                stat="density",
                element="step",
                label=group,
            )

        ax.set_title(f"Distribution of {marker}")
        ax.set_xlabel(f"{marker}")
        ax.legend()

        # Adjust x-axis limits
        side_offset = 0.05 * (upper - lower)
        ax.set_xlim(lower - side_offset, upper + side_offset)

    # Plot each marker
    plot_histogram(axs[0], "Log Free β-hCG")
    plot_histogram(axs[1], "Log PAPP-A")
    plot_histogram(axs[2], "Log NT")

    plt.tight_layout()
    plt.show()

    # Print summary statistics
    print("\nSummary Statistics (Log MoM):")
    for marker in ["Log Free β-hCG", "Log PAPP-A", "Log NT"]:
        print(f"\n{marker}:")
        for group in ["Down Syndrome", "Control"]:
            data = df[df["Group"] == group][marker]
            print(
                f"  {group} - Median: {data.median():.2f}, "
                f"Mean: {data.mean():.2f}, "
                f"SD: {data.std():.2f}"
            )

    # Print summary statistics in original scale (MoM)
    print("\nSummary Statistics (MoM):")
    for marker in ["Log Free β-hCG", "Log PAPP-A", "Log NT"]:
        print(f"\n{marker[4:]}:")  # Remove 'Log ' prefix
        for group in ["Down Syndrome", "Control"]:
            data = np.exp(df[df["Group"] == group][marker])
            print(
                f"  {group} - Median: {data.median():.2f}, "
                f"Mean: {data.mean():.2f}, "
                f"SD: {data.std():.2f}"
            )


if __name__ == "__main__":
    generate_and_visualize_population()
