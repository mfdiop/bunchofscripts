"""
Example Simulations for IBD Method Evaluation

This script demonstrates various simulation scenarios for evaluating
IBD detection methods in P. falciparum.

Author: Generated for IBD method evaluation
Date: 2025-10-01
"""

import numpy as np
from simulation_runner import SimulationRunner


def example_1_basic_simulation():
    """
    Example 1: Basic single population simulation.
    
    Simulates a simple panmictic population to establish baseline IBD patterns.
    """
    print("\n" + "=" * 70)
    print("EXAMPLE 1: Basic Single Population Simulation")
    print("=" * 70)
    
    runner = SimulationRunner(output_dir="results_example1")
    
    # Single population of 100 individuals
    runner.setup_single_population(pop_size=100, recombination_rate=7.4e-7)
    
    # Run for 20 generations
    individuals = runner.run_simulation(num_generations=20, allow_selfing=True)
    
    # Sample 30 individuals for detailed analysis
    sampled = list(np.random.choice(individuals, size=30, replace=False))
    
    # Compute ground truth IBD (only segments > 100kb)
    true_ibd = runner.compute_ground_truth_ibd(sampled, min_length=100000)
    
    # Generate genotype data with 1000 SNPs per chromosome
    genotype_data = runner.generate_genotype_data(sampled, snps_per_chrom=1000)
    
    # Call IBD using all methods
    called_ibd = runner.call_ibd_all_methods(genotype_data, sampled)
    
    # Evaluate methods
    evaluation_df = runner.evaluate_methods()
    
    # Save results
    runner.save_results(prefix="basic_simulation")
    ibd_df = runner.export_ibd_for_plotting(sampled, 
                                            "results_example1/ibd_proportions.csv")
    
    print("\nResults saved to results_example1/")
    print(f"Average IBD proportion: {ibd_df['ibd_proportion'].mean():.4f}")
    
    return runner


def example_2_with_bottleneck():
    """
    Example 2: Population with bottleneck event.
    
    Simulates a bottleneck to increase IBD sharing, mimicking epidemic scenarios.
    """
    print("\n" + "=" * 70)
    print("EXAMPLE 2: Population with Bottleneck")
    print("=" * 70)
    
    runner = SimulationRunner(output_dir="results_example2")
    
    # Start with population of 200
    runner.setup_single_population(pop_size=200, recombination_rate=7.4e-7)
    
    # Define bottleneck at generation 10 (reduce to 20 individuals)
    demographic_events = [
        (10, 'bottleneck', {'population': 0, 'new_size': 20}),
        (11, 'size_change', {'population': 0, 'new_size': 100})  # Expand after bottleneck
    ]
    
    # Run for 30 generations with bottleneck
    individuals = runner.run_simulation(
        num_generations=30,
        allow_selfing=True,
        demographic_events=demographic_events
    )
    
    # Sample individuals
    sampled = list(np.random.choice(individuals, size=30, replace=False))
    
    # Compute ground truth IBD
    true_ibd = runner.compute_ground_truth_ibd(sampled, min_length=100000)
    
    # Generate genotype data
    genotype_data = runner.generate_genotype_data(sampled, snps_per_chrom=1000)
    
    # Call IBD using all methods
    called_ibd = runner.call_ibd_all_methods(genotype_data, sampled)
    
    # Evaluate methods
    evaluation_df = runner.evaluate_methods()
    
    # Save results
    runner.save_results(prefix="bottleneck_simulation")
    ibd_df = runner.export_ibd_for_plotting(sampled,
                                            "results_example2/ibd_proportions.csv")
    
    print("\nResults saved to results_example2/")
    print(f"Average IBD proportion: {ibd_df['ibd_proportion'].mean():.4f}")
    print("(Expected to be higher due to bottleneck)")
    
    return runner


def example_3_two_populations_with_migration():
    """
    Example 3: Two populations with migration.
    
    Simulates population structure with varying levels of gene flow.
    """
    print("\n" + "=" * 70)
    print("EXAMPLE 3: Two Populations with Migration")
    print("=" * 70)
    
    runner = SimulationRunner(output_dir="results_example3")
    
    # Two populations: 80 and 60 individuals
    pop_sizes = [80, 60]
    
    # Migration matrix: 5% migration per generation
    migration_rate = 0.05
    migration_matrix = np.array([
        [1 - migration_rate, migration_rate],      # Pop 0
        [migration_rate, 1 - migration_rate]       # Pop 1
    ])
    
    runner.setup_multi_population(
        pop_sizes=pop_sizes,
        migration_matrix=migration_matrix,
        recombination_rate=7.4e-7
    )
    
    # Run for 25 generations
    individuals = runner.run_simulation(num_generations=25, allow_selfing=True)
    
    # Sample 20 individuals from each population
    pop0_inds = [ind for ind in individuals if ind.population == 0]
    pop1_inds = [ind for ind in individuals if ind.population == 1]
    
    sampled_pop0 = list(np.random.choice(pop0_inds, size=min(20, len(pop0_inds)), replace=False))
    sampled_pop1 = list(np.random.choice(pop1_inds, size=min(20, len(pop1_inds)), replace=False))
    sampled = sampled_pop0 + sampled_pop1
    
    # Compute ground truth IBD
    true_ibd = runner.compute_ground_truth_ibd(sampled, min_length=100000)
    
    # Generate genotype data
    genotype_data = runner.generate_genotype_data(sampled, snps_per_chrom=1000)
    
    # Call IBD using all methods
    called_ibd = runner.call_ibd_all_methods(genotype_data, sampled)
    
    # Evaluate methods
    evaluation_df = runner.evaluate_methods()
    
    # Save results
    runner.save_results(prefix="migration_simulation")
    ibd_df = runner.export_ibd_for_plotting(sampled,
                                            "results_example3/ibd_proportions.csv")
    
    # Analyze within vs between population IBD
    within_pop = ibd_df[ibd_df['same_population'] == True]
    between_pop = ibd_df[ibd_df['same_population'] == False]
    
    print("\nResults saved to results_example3/")
    print(f"Within-population IBD: {within_pop['ibd_proportion'].mean():.4f}")
    print(f"Between-population IBD: {between_pop['ibd_proportion'].mean():.4f}")
    print("(Within-population IBD should be higher)")
    
    return runner


def example_4_isolated_populations():
    """
    Example 4: Two isolated populations (no migration).
    
    Creates distinct populations to test method performance with population structure.
    """
    print("\n" + "=" * 70)
    print("EXAMPLE 4: Two Isolated Populations")
    print("=" * 70)
    
    runner = SimulationRunner(output_dir="results_example4")
    
    # Two populations: 50 each, no migration
    pop_sizes = [50, 50]
    
    # No migration (diagonal matrix)
    migration_matrix = np.array([
        [1.0, 0.0],
        [0.0, 1.0]
    ])
    
    runner.setup_multi_population(
        pop_sizes=pop_sizes,
        migration_matrix=migration_matrix,
        recombination_rate=7.4e-7
    )
    
    # Run for 15 generations
    individuals = runner.run_simulation(num_generations=15, allow_selfing=True)
    
    # Sample from both populations
    pop0_inds = [ind for ind in individuals if ind.population == 0]
    pop1_inds = [ind for ind in individuals if ind.population == 1]
    
    sampled_pop0 = list(np.random.choice(pop0_inds, size=min(15, len(pop0_inds)), replace=False))
    sampled_pop1 = list(np.random.choice(pop1_inds, size=min(15, len(pop1_inds)), replace=False))
    sampled = sampled_pop0 + sampled_pop1
    
    # Compute ground truth IBD
    true_ibd = runner.compute_ground_truth_ibd(sampled, min_length=100000)
    
    # Generate genotype data
    genotype_data = runner.generate_genotype_data(sampled, snps_per_chrom=1000)
    
    # Call IBD using all methods
    called_ibd = runner.call_ibd_all_methods(genotype_data, sampled)
    
    # Evaluate methods
    evaluation_df = runner.evaluate_methods()
    
    # Save results
    runner.save_results(prefix="isolated_simulation")
    ibd_df = runner.export_ibd_for_plotting(sampled,
                                            "results_example4/ibd_proportions.csv")
    
    # Analyze within vs between population IBD
    within_pop = ibd_df[ibd_df['same_population'] == True]
    between_pop = ibd_df[ibd_df['same_population'] == False]
    
    print("\nResults saved to results_example4/")
    print(f"Within-population IBD: {within_pop['ibd_proportion'].mean():.4f}")
    print(f"Between-population IBD: {between_pop['ibd_proportion'].mean():.4f}")
    print("(Between-population IBD should be near zero)")
    
    return runner


def example_5_varying_snp_density():
    """
    Example 5: Compare method performance with varying SNP density.
    
    Tests how SNP density affects IBD calling accuracy.
    """
    print("\n" + "=" * 70)
    print("EXAMPLE 5: Varying SNP Density")
    print("=" * 70)
    
    snp_densities = [100, 500, 1000, 2000]
    results_summary = []
    
    for snp_count in snp_densities:
        print(f"\n--- Testing with {snp_count} SNPs per chromosome ---")
        
        runner = SimulationRunner(output_dir=f"results_example5_snps{snp_count}")
        
        # Single population
        runner.setup_single_population(pop_size=80, recombination_rate=7.4e-7)
        
        # Run for 15 generations
        individuals = runner.run_simulation(num_generations=15, allow_selfing=True)
        
        # Sample individuals
        sampled = list(np.random.choice(individuals, size=25, replace=False))
        
        # Compute ground truth IBD
        true_ibd = runner.compute_ground_truth_ibd(sampled, min_length=100000)
        
        # Generate genotype data with varying SNP density
        genotype_data = runner.generate_genotype_data(sampled, snps_per_chrom=snp_count)
        
        # Call IBD using all methods
        called_ibd = runner.call_ibd_all_methods(genotype_data, sampled)
        
        # Evaluate methods
        evaluation_df = runner.evaluate_methods()
        
        # Save results
        runner.save_results(prefix=f"snp_density_{snp_count}")
        
        # Compute summary
        summary = evaluation_df.groupby('method')[['sensitivity', 'precision', 'f1_score']].mean()
        summary['snp_count'] = snp_count
        results_summary.append(summary)
        
        print(f"Mean F1 scores:")
        for method in summary.index:
            print(f"  {method}: {summary.loc[method, 'f1_score']:.4f}")
    
    print("\n" + "=" * 70)
    print("SNP Density Comparison Complete")
    print("=" * 70)
    
    return results_summary


def run_all_examples():
    """
    Run all example simulations.
    """
    print("\n" + "#" * 70)
    print("# Running All IBD Simulation Examples")
    print("#" * 70)
    
    examples = [
        ("Basic Simulation", example_1_basic_simulation),
        ("Bottleneck", example_2_with_bottleneck),
        ("Migration", example_3_two_populations_with_migration),
        ("Isolated Populations", example_4_isolated_populations),
        ("SNP Density Comparison", example_5_varying_snp_density),
    ]
    
    results = {}
    
    for name, example_func in examples:
        try:
            print(f"\n\nRunning: {name}")
            result = example_func()
            results[name] = result
            print(f"✓ {name} completed successfully")
        except Exception as e:
            print(f"✗ {name} failed with error: {e}")
            import traceback
            traceback.print_exc()
    
    print("\n" + "#" * 70)
    print("# All Examples Complete!")
    print("#" * 70)
    print("\nResults are saved in respective directories:")
    print("  - results_example1/ : Basic simulation")
    print("  - results_example2/ : Bottleneck simulation")
    print("  - results_example3/ : Migration simulation")
    print("  - results_example4/ : Isolated populations")
    print("  - results_example5_snps*/ : SNP density comparison")
    
    return results


if __name__ == "__main__":
    # Run all examples
    results = run_all_examples()
    
    # Or run individual examples:
    # example_1_basic_simulation()
    # example_2_with_bottleneck()
    # example_3_two_populations_with_migration()
    # example_4_isolated_populations()
    # example_5_varying_snp_density()
