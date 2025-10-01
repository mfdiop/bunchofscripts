"""
Simulation Runner and Output Tools

Provides high-level functions to run simulations and export results for analysis.

Author: Generated for IBD method evaluation
Date: 2025-10-01
"""

import numpy as np
import pandas as pd
from typing import List, Dict, Tuple, Optional
import json
from pathlib import Path
import pickle

from ibd_simulation import (
    Simulator, ChromosomeStructure, Individual,
    summarize_ibd
)
from ibd_calling import (
    ThresholdIBDCaller, HMMIBDCaller, WindowIBDCaller,
    CalledIBDSegment, evaluate_ibd_calling
)


class SimulationRunner:
    """
    High-level interface for running IBD simulations and evaluations.
    """
    
    def __init__(self, output_dir: str = "simulation_results"):
        """
        Initialize simulation runner.
        
        Args:
            output_dir: Directory to save results
        """
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True, parents=True)
        self.simulator = None
        self.results = {}
    
    def setup_single_population(self, pop_size: int = 100,
                               recombination_rate: float = 7.4e-7) -> Simulator:
        """
        Set up a single population simulation.
        
        Args:
            pop_size: Population size
            recombination_rate: Recombination rate per bp per generation
            
        Returns:
            Configured simulator
        """
        chrom_structure = ChromosomeStructure(recombination_rate)
        self.simulator = Simulator(chrom_structure)
        self.simulator.add_population(0, pop_size)
        return self.simulator
    
    def setup_multi_population(self, pop_sizes: List[int],
                              migration_matrix: Optional[np.ndarray] = None,
                              recombination_rate: float = 7.4e-7) -> Simulator:
        """
        Set up a multi-population simulation.
        
        Args:
            pop_sizes: List of population sizes
            migration_matrix: Migration matrix (default: no migration)
            recombination_rate: Recombination rate per bp per generation
            
        Returns:
            Configured simulator
        """
        chrom_structure = ChromosomeStructure(recombination_rate)
        self.simulator = Simulator(chrom_structure)
        
        for pop_id, size in enumerate(pop_sizes):
            self.simulator.add_population(pop_id, size)
        
        if migration_matrix is not None:
            self.simulator.set_migration_matrix(migration_matrix)
        
        return self.simulator
    
    def run_simulation(self, num_generations: int,
                      allow_selfing: bool = False,
                      demographic_events: Optional[List[Tuple[int, str, Dict]]] = None,
                      sample_size: Optional[int] = None) -> List[Individual]:
        """
        Run the simulation and return sampled individuals.
        
        Args:
            num_generations: Number of generations to simulate
            allow_selfing: Whether to allow self-fertilization
            demographic_events: List of demographic events
            sample_size: Number of individuals to sample (None = all)
            
        Returns:
            List of sampled individuals
        """
        if self.simulator is None:
            raise ValueError("Must call setup_single_population or setup_multi_population first")
        
        print(f"Running simulation for {num_generations} generations...")
        self.simulator.run_simulation(num_generations, allow_selfing, demographic_events)
        
        # Get all individuals
        all_individuals = self.simulator.get_all_individuals()
        
        # Sample if requested
        if sample_size is not None and sample_size < len(all_individuals):
            sampled = list(np.random.choice(all_individuals, size=sample_size, replace=False))
        else:
            sampled = all_individuals
        
        print(f"Simulation complete. Sampled {len(sampled)} individuals.")
        
        self.results['individuals'] = sampled
        return sampled
    
    def compute_ground_truth_ibd(self, individuals: Optional[List[Individual]] = None,
                                min_length: float = 0) -> Dict[Tuple[int, int], List[Tuple[int, float, float]]]:
        """
        Compute ground truth IBD for sampled individuals.
        
        Args:
            individuals: List of individuals (default: from run_simulation)
            min_length: Minimum IBD segment length to report
            
        Returns:
            Dictionary of pairwise IBD segments
        """
        if individuals is None:
            individuals = self.results.get('individuals')
            if individuals is None:
                raise ValueError("Must run simulation first or provide individuals")
        
        print(f"Computing ground truth IBD for {len(individuals)} individuals...")
        ibd_results = self.simulator.compute_pairwise_ibd(individuals, min_length)
        
        self.results['true_ibd'] = ibd_results
        print(f"Found IBD in {len(ibd_results)} pairs.")
        
        return ibd_results
    
    def generate_genotype_data(self, individuals: List[Individual],
                              snps_per_chrom: int = 1000,
                              mutation_rate: float = 1e-8) -> Dict[int, np.ndarray]:
        """
        Generate SNP genotype data from individuals.
        
        This simulates genotyping by placing SNPs and determining genotypes
        based on the underlying IBD structure.
        
        Args:
            individuals: List of individuals
            snps_per_chrom: Number of SNPs per chromosome
            mutation_rate: Mutation rate (for creating variation)
            
        Returns:
            Dictionary mapping chromosome index to genotype matrix
        """
        print(f"Generating genotype data with {snps_per_chrom} SNPs per chromosome...")
        
        genotype_data = {}
        num_chroms = self.simulator.chrom_structure.num_chromosomes
        
        for chrom_idx in range(num_chroms):
            chrom_length = self.simulator.chrom_structure.lengths[chrom_idx]
            
            # Place SNPs uniformly
            snp_positions = np.linspace(0, chrom_length, snps_per_chrom)
            
            # Initialize genotype matrix
            # For each individual, we'll use the first haplotype for simplicity
            # (or we could encode both haplotypes)
            genotypes = np.zeros((snps_per_chrom, len(individuals)), dtype=int)
            
            # Assign genotypes based on founder haplotypes
            # First, assign random alleles to each founder at each SNP
            max_founder = 0
            for ind in individuals:
                for chrom_pair in ind.chromosomes:
                    for chrom in chrom_pair:
                        for seg in chrom.segments:
                            max_founder = max(max_founder, seg.founder_id)
            
            # Create founder alleles at each SNP with mutation
            founder_alleles = {}
            for founder_id in range(max_founder + 1):
                # Random founder alleles with given allele frequency
                base_allele = np.random.binomial(1, 0.5)
                # Add some variation across SNPs
                alleles = np.random.binomial(1, 0.1, snps_per_chrom)
                founder_alleles[founder_id] = (alleles == 1).astype(int)
                # Set a baseline
                founder_alleles[founder_id][0] = base_allele
            
            # Assign genotypes to individuals based on their segments
            for ind_idx, ind in enumerate(individuals):
                maternal_chrom, paternal_chrom = ind.chromosomes[chrom_idx]
                
                # Use maternal haplotype for genotype
                for snp_idx, snp_pos in enumerate(snp_positions):
                    # Find which segment this SNP falls in
                    for seg in maternal_chrom.segments:
                        if seg.start <= snp_pos < seg.end:
                            genotypes[snp_idx, ind_idx] = founder_alleles[seg.founder_id][snp_idx]
                            break
            
            genotype_data[chrom_idx] = genotypes
        
        self.results['genotype_data'] = genotype_data
        print("Genotype data generated.")
        
        return genotype_data
    
    def call_ibd_all_methods(self, genotype_data: Dict[int, np.ndarray],
                            individuals: List[Individual]) -> Dict[str, Dict]:
        """
        Call IBD using all implemented methods.
        
        Args:
            genotype_data: Genotype data from generate_genotype_data
            individuals: List of individuals
            
        Returns:
            Dictionary mapping method name to results
        """
        print("Calling IBD using multiple methods...")
        
        chrom_lengths = self.simulator.chrom_structure.lengths
        
        # Initialize callers
        callers = {
            'threshold': ThresholdIBDCaller(),
            'hmm': HMMIBDCaller(),
            'window': WindowIBDCaller(),
        }
        
        all_results = {}
        
        for method_name, caller in callers.items():
            print(f"  Running {method_name} method...")
            method_results = {}
            
            # Call IBD for all pairs
            for i in range(len(individuals)):
                for j in range(i + 1, len(individuals)):
                    called_segs = caller.call_ibd(genotype_data, i, j, chrom_lengths)
                    
                    if called_segs:
                        ind1_id = individuals[i].id
                        ind2_id = individuals[j].id
                        method_results[(ind1_id, ind2_id)] = called_segs
            
            all_results[method_name] = method_results
            print(f"    Found IBD in {len(method_results)} pairs")
        
        self.results['called_ibd'] = all_results
        
        return all_results
    
    def evaluate_methods(self, true_ibd: Optional[Dict] = None,
                        called_ibd: Optional[Dict] = None) -> pd.DataFrame:
        """
        Evaluate all IBD calling methods against ground truth.
        
        Args:
            true_ibd: Ground truth IBD (default: from results)
            called_ibd: Called IBD (default: from results)
            
        Returns:
            DataFrame with evaluation metrics for each method
        """
        if true_ibd is None:
            true_ibd = self.results.get('true_ibd')
        if called_ibd is None:
            called_ibd = self.results.get('called_ibd')
        
        if true_ibd is None or called_ibd is None:
            raise ValueError("Must compute ground truth and call IBD first")
        
        print("Evaluating IBD calling methods...")
        
        total_genome_length = sum(self.simulator.chrom_structure.lengths)
        
        eval_results = []
        
        for method_name, method_results in called_ibd.items():
            print(f"  Evaluating {method_name}...")
            
            method_metrics = []
            
            # Get all pairs that have either true or called IBD
            all_pairs = set(true_ibd.keys()) | set(method_results.keys())
            
            for pair in all_pairs:
                true_segs = true_ibd.get(pair, [])
                called_segs = method_results.get(pair, [])
                
                # Convert called segments to same format as true
                called_segs_formatted = [
                    (seg.chrom_idx, seg.start, seg.end) for seg in called_segs
                ]
                
                metrics = evaluate_ibd_calling(
                    true_segs, called_segs, total_genome_length
                )
                metrics['method'] = method_name
                metrics['pair'] = pair
                method_metrics.append(metrics)
            
            eval_results.extend(method_metrics)
        
        df = pd.DataFrame(eval_results)
        
        # Compute summary statistics
        summary = df.groupby('method')[['sensitivity', 'precision', 'f1_score']].agg(['mean', 'std'])
        print("\nSummary by method:")
        print(summary)
        
        self.results['evaluation'] = df
        
        return df
    
    def save_results(self, prefix: str = "simulation"):
        """
        Save all results to disk.
        
        Args:
            prefix: Prefix for output files
        """
        print(f"Saving results to {self.output_dir}...")
        
        # Save individuals (pickle for full object)
        if 'individuals' in self.results:
            with open(self.output_dir / f"{prefix}_individuals.pkl", 'wb') as f:
                pickle.dump(self.results['individuals'], f)
        
        # Save true IBD (JSON format)
        if 'true_ibd' in self.results:
            true_ibd_serializable = {
                f"{k[0]}_{k[1]}": v for k, v in self.results['true_ibd'].items()
            }
            with open(self.output_dir / f"{prefix}_true_ibd.json", 'w') as f:
                json.dump(true_ibd_serializable, f, indent=2)
        
        # Save genotype data
        if 'genotype_data' in self.results:
            for chrom_idx, genotypes in self.results['genotype_data'].items():
                np.save(
                    self.output_dir / f"{prefix}_genotypes_chr{chrom_idx}.npy",
                    genotypes
                )
        
        # Save evaluation results
        if 'evaluation' in self.results:
            self.results['evaluation'].to_csv(
                self.output_dir / f"{prefix}_evaluation.csv",
                index=False
            )
        
        # Save summary statistics
        if 'evaluation' in self.results:
            summary = self.results['evaluation'].groupby('method')[
                ['sensitivity', 'precision', 'f1_score']
            ].agg(['mean', 'std', 'median'])
            summary.to_csv(self.output_dir / f"{prefix}_summary.csv")
        
        print("Results saved successfully.")
    
    def export_ibd_for_plotting(self, individuals: Optional[List[Individual]] = None,
                               output_file: Optional[str] = None) -> pd.DataFrame:
        """
        Export IBD data in format suitable for plotting.
        
        Args:
            individuals: List of individuals (default: from results)
            output_file: Optional output CSV file
            
        Returns:
            DataFrame with IBD proportions
        """
        if individuals is None:
            individuals = self.results.get('individuals')
        
        true_ibd = self.results.get('true_ibd', {})
        
        rows = []
        
        for i, ind1 in enumerate(individuals):
            for j, ind2 in enumerate(individuals):
                if i >= j:
                    continue
                
                pair = (ind1.id, ind2.id)
                segs = true_ibd.get(pair, [])
                
                total_ibd = sum(end - start for _, start, end in segs)
                total_genome = sum(self.simulator.chrom_structure.lengths)
                ibd_proportion = total_ibd / total_genome
                
                rows.append({
                    'ind1_id': ind1.id,
                    'ind2_id': ind2.id,
                    'ind1_pop': ind1.population,
                    'ind2_pop': ind2.population,
                    'same_population': ind1.population == ind2.population,
                    'ibd_proportion': ibd_proportion,
                    'ibd_length_mb': total_ibd / 1e6,
                    'num_segments': len(segs)
                })
        
        df = pd.DataFrame(rows)
        
        if output_file:
            df.to_csv(output_file, index=False)
        
        return df


def run_simple_example():
    """
    Run a simple example simulation.
    """
    print("=" * 60)
    print("Running simple IBD simulation example")
    print("=" * 60)
    
    runner = SimulationRunner(output_dir="example_results")
    
    # Set up single population
    runner.setup_single_population(pop_size=50, recombination_rate=7.4e-7)
    
    # Run for 10 generations
    individuals = runner.run_simulation(num_generations=10, allow_selfing=True)
    
    # Sample subset for analysis
    sampled = list(np.random.choice(individuals, size=min(20, len(individuals)), replace=False))
    
    # Compute ground truth IBD
    true_ibd = runner.compute_ground_truth_ibd(sampled, min_length=50000)
    
    # Generate genotype data
    genotype_data = runner.generate_genotype_data(sampled, snps_per_chrom=500)
    
    # Call IBD using all methods
    called_ibd = runner.call_ibd_all_methods(genotype_data, sampled)
    
    # Evaluate methods
    evaluation_df = runner.evaluate_methods()
    
    # Save results
    runner.save_results(prefix="simple_example")
    
    # Export for plotting
    ibd_df = runner.export_ibd_for_plotting(sampled, "example_results/ibd_proportions.csv")
    
    print("\n" + "=" * 60)
    print("Example complete! Results saved to example_results/")
    print("=" * 60)
    
    return runner


if __name__ == "__main__":
    run_simple_example()
