import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import pearsonr
import sys

def load_data(filename):
    data = np.loadtxt(filename)
    costs = data[:, 0]
    sim_edge_best = data[:, 1]
    sim_cycle_best = data[:, 2]
    avg_sim_edge_all = data[:, 3]
    avg_sim_cycle_all = data[:, 4]
    return costs, sim_edge_best, sim_cycle_best, avg_sim_edge_all, avg_sim_cycle_all

def plot_and_corr(x, y, title, xlabel, ylabel, subplot_position):
    plt.subplot(2, 2, subplot_position)

    plt.scatter(x, y, s=10, alpha=0.6)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    
    corr, _ = pearsonr(x, y)
    plt.title(f"Korelacja = {corr:.4f}")
    plt.grid(True)

def main(filename):
    costs, sim_edge_best, sim_cycle_best, avg_sim_edge_all, avg_sim_cycle_all = load_data(filename)

    plt.figure(figsize=(14, 10))

    plot_and_corr(costs, sim_edge_best, "krawędź; Best",
                  "Funkcja celu", "Wspólne krawędzie", 1)

    plot_and_corr(costs, sim_cycle_best, "cykl: Best",
                  "Funkcja celu", "Wspólne pary wierzchołków", 2)

    plot_and_corr(costs, avg_sim_edge_all, "krawędź; Mean",
                  "Funkcja celu", "Średnie wspólne krawędzie", 3)

    plot_and_corr(costs, avg_sim_cycle_all, "cykl: Mean",
                  "Funkcja celu", "Średnie wspólne pary wierzchołków", 4)

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Użycie: python analyze_plot.py <plik_danych.txt>")
        sys.exit(1)

    main(sys.argv[1])
