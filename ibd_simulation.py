"""
IBD Segment Tracking Simulator for P. falciparum

This module implements a forward-time population genetics simulator that explicitly
tracks Identity-by-Descent (IBD) segments through recombination and inheritance.

Key Features:
- Proper chromosome structure with recombination
- Explicit IBD segment tracking with founder labels
- Multiple populations with migration
- Demographic events (bottlenecks, population size changes)
- Ground truth IBD for method evaluation

Author: Generated for IBD method evaluation
Date: 2025-10-01
"""

import numpy as np
from dataclasses import dataclass, field
from typing import List, Dict, Tuple, Optional, Set
from collections import defaultdict
import warnings


@dataclass
class Segment:
    """
    Represents a chromosomal segment with IBD information.
    
    Attributes:
        start: Start position in base pairs
        end: End position in base pairs
        founder_id: Original founder haplotype identifier
        generation_origin: Generation when this segment was created
    """
    start: float
    end: float
    founder_id: int
    generation_origin: int = 0
    
    def length(self) -> float:
        """Return the length of this segment in base pairs."""
        return self.end - self.start
    
    def overlaps(self, other: 'Segment') -> bool:
        """Check if this segment overlaps with another."""
        return not (self.end <= other.start or self.start >= other.end)
    
    def __repr__(self) -> str:
        return f"Segment({self.start:.0f}-{self.end:.0f}, founder={self.founder_id})"


@dataclass
class Chromosome:
    """
    Represents a haploid chromosome as a list of IBD segments.
    
    Attributes:
        segments: List of non-overlapping segments covering the chromosome
        length: Total chromosome length in base pairs
    """
    segments: List[Segment]
    length: float
    
    def __post_init__(self):
        """Validate that segments are non-overlapping and sorted."""
        if not self.segments:
            raise ValueError("Chromosome must have at least one segment")
        
        # Sort segments by start position
        self.segments.sort(key=lambda s: s.start)
        
        # Verify no gaps or overlaps
        for i in range(len(self.segments) - 1):
            if self.segments[i].end != self.segments[i+1].start:
                raise ValueError(f"Segments must be contiguous: gap between "
                               f"{self.segments[i]} and {self.segments[i+1]}")
    
    def copy(self) -> 'Chromosome':
        """Create a deep copy of this chromosome."""
        return Chromosome(
            segments=[Segment(s.start, s.end, s.founder_id, s.generation_origin) 
                     for s in self.segments],
            length=self.length
        )
    
    def get_ibd_with(self, other: 'Chromosome', min_length: float = 0) -> List[Tuple[float, float]]:
        """
        Find IBD segments shared with another chromosome.
        
        Args:
            other: Another chromosome to compare against
            min_length: Minimum segment length to report
            
        Returns:
            List of (start, end) tuples representing shared IBD regions
        """
        ibd_regions = []
        i, j = 0, 0
        
        while i < len(self.segments) and j < len(other.segments):
            seg1 = self.segments[i]
            seg2 = other.segments[j]
            
            # Check if same founder
            if seg1.founder_id == seg2.founder_id:
                # Find overlap
                start = max(seg1.start, seg2.start)
                end = min(seg1.end, seg2.end)
                
                if end > start and (end - start) >= min_length:
                    ibd_regions.append((start, end))
            
            # Move to next segment
            if seg1.end < seg2.end:
                i += 1
            elif seg2.end < seg1.end:
                j += 1
            else:
                i += 1
                j += 1
        
        return ibd_regions


@dataclass
class Individual:
    """
    Represents a diploid individual with multiple chromosomes.
    
    Attributes:
        id: Unique identifier
        chromosomes: List of [maternal, paternal] chromosome pairs
        generation: Generation number
        population: Population identifier
    """
    id: int
    chromosomes: List[Tuple[Chromosome, Chromosome]]
    generation: int
    population: int = 0
    
    def get_haplotype(self, chrom_idx: int, parental: int) -> Chromosome:
        """
        Get a specific haplotype.
        
        Args:
            chrom_idx: Chromosome index
            parental: 0 for maternal, 1 for paternal
            
        Returns:
            The requested chromosome
        """
        return self.chromosomes[chrom_idx][parental]


class ChromosomeStructure:
    """
    Defines the chromosome structure for P. falciparum.
    
    P. falciparum has 14 chromosomes with varying lengths and recombination rates.
    """
    
    # P. falciparum 3D7 chromosome lengths (approximate, in base pairs)
    CHROMOSOME_LENGTHS = [
        640851,   # Chr 1
        947102,   # Chr 2
        1067971,  # Chr 3
        1200490,  # Chr 4
        1343557,  # Chr 5
        1418242,  # Chr 6
        1445207,  # Chr 7
        1472805,  # Chr 8
        1541735,  # Chr 9
        1687656,  # Chr 10
        2038340,  # Chr 11
        2271494,  # Chr 12
        2925236,  # Chr 13
        3291936,  # Chr 14
    ]
    
    def __init__(self, recombination_rate: float = 7.4e-7):
        """
        Initialize chromosome structure.
        
        Args:
            recombination_rate: Recombination rate per base pair per generation
                              (default ~7.4e-7 based on P. falciparum estimates)
        """
        self.lengths = self.CHROMOSOME_LENGTHS
        self.num_chromosomes = len(self.lengths)
        self.recombination_rate = recombination_rate
    
    def get_recombination_breakpoints(self, chrom_idx: int) -> List[float]:
        """
        Generate recombination breakpoints for a chromosome.
        
        Uses a Poisson process to determine the number of crossovers,
        then places them uniformly along the chromosome.
        
        Args:
            chrom_idx: Index of the chromosome
            
        Returns:
            Sorted list of breakpoint positions
        """
        length = self.lengths[chrom_idx]
        
        # Expected number of crossovers
        expected_crossovers = length * self.recombination_rate
        
        # Sample number of crossovers from Poisson distribution
        num_crossovers = np.random.poisson(expected_crossovers)
        
        if num_crossovers == 0:
            return []
        
        # Sample crossover positions uniformly
        breakpoints = sorted(np.random.uniform(0, length, num_crossovers))
        
        return breakpoints


class Population:
    """
    Manages a population of individuals with demographic events.
    
    Attributes:
        individuals: Current individuals in the population
        size: Current population size
        id: Population identifier
    """
    
    def __init__(self, pop_id: int, initial_size: int, 
                 chrom_structure: ChromosomeStructure):
        """
        Initialize a population.
        
        Args:
            pop_id: Population identifier
            initial_size: Initial population size
            chrom_structure: Chromosome structure definition
        """
        self.id = pop_id
        self.size = initial_size
        self.chrom_structure = chrom_structure
        self.individuals: List[Individual] = []
        self.generation = 0
        self._next_individual_id = 0
        self._next_founder_id = 0
    
    def initialize_founders(self):
        """
        Create founder individuals with unique founder labels.
        
        Each founder haplotype gets a unique label for IBD tracking.
        """
        for i in range(self.size):
            chromosomes = []
            
            for chrom_idx in range(self.chrom_structure.num_chromosomes):
                length = self.chrom_structure.lengths[chrom_idx]
                
                # Create maternal chromosome with unique founder
                maternal = Chromosome(
                    segments=[Segment(0, length, self._next_founder_id, 0)],
                    length=length
                )
                self._next_founder_id += 1
                
                # Create paternal chromosome with unique founder
                paternal = Chromosome(
                    segments=[Segment(0, length, self._next_founder_id, 0)],
                    length=length
                )
                self._next_founder_id += 1
                
                chromosomes.append((maternal, paternal))
            
            individual = Individual(
                id=self._next_individual_id,
                chromosomes=chromosomes,
                generation=0,
                population=self.id
            )
            self._next_individual_id += 1
            self.individuals.append(individual)
    
    def recombine(self, parent1_chrom: Chromosome, parent2_chrom: Chromosome,
                  breakpoints: List[float]) -> Chromosome:
        """
        Perform recombination between two parental chromosomes.
        
        Args:
            parent1_chrom: First parental chromosome
            parent2_chrom: Second parental chromosome
            breakpoints: List of recombination breakpoint positions
            
        Returns:
            New recombinant chromosome
        """
        if not breakpoints:
            # No recombination, return copy of one parent
            return parent1_chrom.copy()
        
        new_segments = []
        current_parent = 0  # Start with parent1
        parents = [parent1_chrom, parent2_chrom]
        
        # Add chromosome start
        positions = [0] + sorted(breakpoints) + [parent1_chrom.length]
        
        for i in range(len(positions) - 1):
            start = positions[i]
            end = positions[i + 1]
            
            # Get segments from current parent in this region
            current_chrom = parents[current_parent]
            
            for seg in current_chrom.segments:
                # Find overlap with current region
                if seg.end <= start:
                    continue
                if seg.start >= end:
                    break
                
                # Create new segment for overlapping region
                seg_start = max(seg.start, start)
                seg_end = min(seg.end, end)
                
                new_segments.append(Segment(
                    seg_start, seg_end, seg.founder_id, 
                    self.generation + 1
                ))
            
            # Switch parent at breakpoint
            current_parent = 1 - current_parent
        
        # Merge adjacent segments with same founder
        merged_segments = []
        for seg in new_segments:
            if (merged_segments and 
                merged_segments[-1].founder_id == seg.founder_id and
                merged_segments[-1].end == seg.start):
                # Extend previous segment
                merged_segments[-1].end = seg.end
            else:
                merged_segments.append(seg)
        
        return Chromosome(merged_segments, parent1_chrom.length)
    
    def mate(self, parent1: Individual, parent2: Individual) -> Individual:
        """
        Create offspring from two parents through meiosis and recombination.
        
        Args:
            parent1: First parent
            parent2: Second parent
            
        Returns:
            New offspring individual
        """
        offspring_chromosomes = []
        
        for chrom_idx in range(self.chrom_structure.num_chromosomes):
            # Get parental chromosomes
            p1_maternal, p1_paternal = parent1.chromosomes[chrom_idx]
            p2_maternal, p2_paternal = parent2.chromosomes[chrom_idx]
            
            # Generate recombination breakpoints for each meiosis
            breakpoints1 = self.chrom_structure.get_recombination_breakpoints(chrom_idx)
            breakpoints2 = self.chrom_structure.get_recombination_breakpoints(chrom_idx)
            
            # Create gametes through recombination
            gamete1 = self.recombine(p1_maternal, p1_paternal, breakpoints1)
            gamete2 = self.recombine(p2_maternal, p2_paternal, breakpoints2)
            
            offspring_chromosomes.append((gamete1, gamete2))
        
        offspring = Individual(
            id=self._next_individual_id,
            chromosomes=offspring_chromosomes,
            generation=self.generation + 1,
            population=self.id
        )
        self._next_individual_id += 1
        
        return offspring
    
    def advance_generation(self, allow_selfing: bool = False):
        """
        Advance the population by one generation using Wright-Fisher model.
        
        Args:
            allow_selfing: Whether to allow self-fertilization
        """
        new_individuals = []
        
        for _ in range(self.size):
            # Sample two parents with replacement
            if allow_selfing:
                parent1 = np.random.choice(self.individuals)
                parent2 = np.random.choice(self.individuals)
            else:
                # Sample without replacement if selfing not allowed
                parents = np.random.choice(self.individuals, size=2, replace=False)
                parent1, parent2 = parents[0], parents[1]
            
            offspring = self.mate(parent1, parent2)
            new_individuals.append(offspring)
        
        self.individuals = new_individuals
        self.generation += 1
    
    def apply_bottleneck(self, new_size: int):
        """
        Apply a population bottleneck by randomly sampling individuals.
        
        Args:
            new_size: New population size after bottleneck
        """
        if new_size > len(self.individuals):
            warnings.warn(f"Bottleneck size {new_size} larger than current size "
                         f"{len(self.individuals)}. No bottleneck applied.")
            return
        
        self.individuals = list(np.random.choice(
            self.individuals, size=new_size, replace=False
        ))
        self.size = new_size


class Simulator:
    """
    Main simulator class for multi-population IBD tracking.
    
    Handles multiple populations with migration and demographic events.
    """
    
    def __init__(self, chrom_structure: Optional[ChromosomeStructure] = None):
        """
        Initialize the simulator.
        
        Args:
            chrom_structure: Chromosome structure (default: P. falciparum)
        """
        self.chrom_structure = chrom_structure or ChromosomeStructure()
        self.populations: Dict[int, Population] = {}
        self.generation = 0
        self.migration_matrix: Optional[np.ndarray] = None
    
    def add_population(self, pop_id: int, size: int) -> Population:
        """
        Add a new population to the simulation.
        
        Args:
            pop_id: Population identifier
            size: Population size
            
        Returns:
            The created population
        """
        pop = Population(pop_id, size, self.chrom_structure)
        pop.initialize_founders()
        self.populations[pop_id] = pop
        return pop
    
    def set_migration_matrix(self, migration_matrix: np.ndarray):
        """
        Set migration matrix for multi-population simulations.
        
        Args:
            migration_matrix: Matrix where element [i,j] is the proportion
                            of population i that comes from population j
                            each generation. Rows should sum to 1.
        """
        # Validate matrix
        num_pops = len(self.populations)
        if migration_matrix.shape != (num_pops, num_pops):
            raise ValueError(f"Migration matrix shape {migration_matrix.shape} "
                           f"does not match number of populations {num_pops}")
        
        row_sums = migration_matrix.sum(axis=1)
        if not np.allclose(row_sums, 1.0):
            raise ValueError(f"Migration matrix rows must sum to 1. Got: {row_sums}")
        
        self.migration_matrix = migration_matrix
    
    def advance_generation_with_migration(self, allow_selfing: bool = False):
        """
        Advance all populations by one generation with migration.
        
        Args:
            allow_selfing: Whether to allow self-fertilization
        """
        if self.migration_matrix is None:
            # No migration, just advance each population
            for pop in self.populations.values():
                pop.advance_generation(allow_selfing)
        else:
            # Apply migration
            pop_ids = sorted(self.populations.keys())
            new_individuals_by_pop = {pop_id: [] for pop_id in pop_ids}
            
            for i, pop_id in enumerate(pop_ids):
                pop = self.populations[pop_id]
                
                for _ in range(pop.size):
                    # Determine source populations for parents based on migration
                    source_pop_idx1 = np.random.choice(
                        len(pop_ids), p=self.migration_matrix[i, :]
                    )
                    source_pop_idx2 = np.random.choice(
                        len(pop_ids), p=self.migration_matrix[i, :]
                    )
                    
                    source_pop1 = self.populations[pop_ids[source_pop_idx1]]
                    source_pop2 = self.populations[pop_ids[source_pop_idx2]]
                    
                    # Sample parents
                    parent1 = np.random.choice(source_pop1.individuals)
                    
                    if not allow_selfing and source_pop_idx1 == source_pop_idx2:
                        # Ensure different parents if same population
                        available_parents = [ind for ind in source_pop2.individuals 
                                           if ind.id != parent1.id]
                        if available_parents:
                            parent2 = np.random.choice(available_parents)
                        else:
                            parent2 = np.random.choice(source_pop2.individuals)
                    else:
                        parent2 = np.random.choice(source_pop2.individuals)
                    
                    # Create offspring (using destination population's mate method)
                    offspring = pop.mate(parent1, parent2)
                    offspring.population = pop_id
                    new_individuals_by_pop[pop_id].append(offspring)
            
            # Update populations
            for pop_id in pop_ids:
                self.populations[pop_id].individuals = new_individuals_by_pop[pop_id]
                self.populations[pop_id].generation += 1
        
        self.generation += 1
    
    def run_simulation(self, num_generations: int, allow_selfing: bool = False,
                      demographic_events: Optional[List[Tuple[int, str, Dict]]] = None):
        """
        Run the simulation for multiple generations.
        
        Args:
            num_generations: Number of generations to simulate
            allow_selfing: Whether to allow self-fertilization
            demographic_events: List of (generation, event_type, params) tuples
                              where event_type can be 'bottleneck', 'size_change', etc.
        """
        demographic_events = demographic_events or []
        event_dict = defaultdict(list)
        for gen, event_type, params in demographic_events:
            event_dict[gen].append((event_type, params))
        
        for gen in range(num_generations):
            # Apply demographic events
            if gen in event_dict:
                for event_type, params in event_dict[gen]:
                    if event_type == 'bottleneck':
                        pop_id = params['population']
                        new_size = params['new_size']
                        self.populations[pop_id].apply_bottleneck(new_size)
                    elif event_type == 'size_change':
                        pop_id = params['population']
                        new_size = params['new_size']
                        self.populations[pop_id].size = new_size
            
            # Advance generation
            self.advance_generation_with_migration(allow_selfing)
    
    def get_all_individuals(self) -> List[Individual]:
        """Get all individuals from all populations."""
        individuals = []
        for pop in self.populations.values():
            individuals.extend(pop.individuals)
        return individuals
    
    def compute_pairwise_ibd(self, individuals: Optional[List[Individual]] = None,
                            min_length: float = 0) -> Dict[Tuple[int, int], List[Tuple[int, float, float]]]:
        """
        Compute ground truth pairwise IBD for all individuals.
        
        Args:
            individuals: List of individuals (default: all from all populations)
            min_length: Minimum IBD segment length to report
            
        Returns:
            Dictionary mapping (ind1_id, ind2_id) to list of (chrom_idx, start, end) tuples
        """
        if individuals is None:
            individuals = self.get_all_individuals()
        
        ibd_results = {}
        
        for i in range(len(individuals)):
            for j in range(i + 1, len(individuals)):
                ind1 = individuals[i]
                ind2 = individuals[j]
                
                ibd_segments = []
                
                # Check each chromosome
                for chrom_idx in range(self.chrom_structure.num_chromosomes):
                    # Compare all four haplotype combinations
                    for h1 in [0, 1]:  # ind1 haplotypes
                        for h2 in [0, 1]:  # ind2 haplotypes
                            chrom1 = ind1.get_haplotype(chrom_idx, h1)
                            chrom2 = ind2.get_haplotype(chrom_idx, h2)
                            
                            # Find IBD regions
                            ibd_regions = chrom1.get_ibd_with(chrom2, min_length)
                            
                            for start, end in ibd_regions:
                                ibd_segments.append((chrom_idx, start, end))
                
                if ibd_segments:
                    ibd_results[(ind1.id, ind2.id)] = ibd_segments
        
        return ibd_results


def summarize_ibd(ibd_results: Dict[Tuple[int, int], List[Tuple[int, float, float]]],
                 chrom_structure: ChromosomeStructure) -> Dict[Tuple[int, int], float]:
    """
    Summarize IBD results as genome-wide IBD proportion.
    
    Args:
        ibd_results: Output from compute_pairwise_ibd
        chrom_structure: Chromosome structure for genome size calculation
        
    Returns:
        Dictionary mapping (ind1_id, ind2_id) to proportion of genome IBD
    """
    total_genome_length = sum(chrom_structure.lengths)
    
    summary = {}
    for (ind1_id, ind2_id), segments in ibd_results.items():
        total_ibd_length = sum(end - start for _, start, end in segments)
        ibd_proportion = total_ibd_length / total_genome_length
        summary[(ind1_id, ind2_id)] = ibd_proportion
    
    return summary
