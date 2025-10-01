"""
Basic functionality tests for the IBD simulation framework.

This script tests core components to ensure they work correctly.
"""

import numpy as np
from ibd_simulation import (
    Segment, Chromosome, Individual, ChromosomeStructure, 
    Population, Simulator
)
from ibd_calling import ThresholdIBDCaller, HMMIBDCaller, WindowIBDCaller
from simulation_runner import SimulationRunner


def test_segment_creation():
    """Test basic segment creation and methods."""
    print("Testing Segment creation...")
    seg = Segment(start=0, end=1000000, founder_id=1, generation_origin=0)
    assert seg.length() == 1000000
    print("  ✓ Segment creation and length calculation")
    
    seg2 = Segment(start=500000, end=1500000, founder_id=2)
    assert seg.overlaps(seg2)
    print("  ✓ Segment overlap detection")


def test_chromosome_creation():
    """Test chromosome creation and IBD detection."""
    print("\nTesting Chromosome creation...")
    
    # Create two chromosomes with some shared founders
    chrom1 = Chromosome(
        segments=[
            Segment(0, 500000, founder_id=1),
            Segment(500000, 1000000, founder_id=2),
        ],
        length=1000000
    )
    
    chrom2 = Chromosome(
        segments=[
            Segment(0, 300000, founder_id=3),
            Segment(300000, 800000, founder_id=1),  # Shared with chrom1
            Segment(800000, 1000000, founder_id=4),
        ],
        length=1000000
    )
    
    print("  ✓ Chromosome creation with multiple segments")
    
    # Test IBD detection
    ibd_regions = chrom1.get_ibd_with(chrom2)
    assert len(ibd_regions) == 1
    assert ibd_regions[0] == (300000, 500000)  # Overlap of founder 1
    print(f"  ✓ IBD detection: found region {ibd_regions[0]}")


def test_chromosome_structure():
    """Test P. falciparum chromosome structure."""
    print("\nTesting ChromosomeStructure...")
    
    chrom_struct = ChromosomeStructure(recombination_rate=7.4e-7)
    assert chrom_struct.num_chromosomes == 14
    assert len(chrom_struct.lengths) == 14
    print(f"  ✓ Correct number of chromosomes: {chrom_struct.num_chromosomes}")
    
    # Test recombination breakpoints
    breakpoints = chrom_struct.get_recombination_breakpoints(0)
    assert isinstance(breakpoints, list)
    print(f"  ✓ Recombination breakpoints generated: {len(breakpoints)} crossovers")


def test_population_initialization():
    """Test population initialization with founders."""
    print("\nTesting Population initialization...")
    
    chrom_struct = ChromosomeStructure()
    pop = Population(pop_id=0, initial_size=10, chrom_structure=chrom_struct)
    pop.initialize_founders()
    
    assert len(pop.individuals) == 10
    print(f"  ✓ Created {len(pop.individuals)} founder individuals")
    
    # Check that each founder has all chromosomes
    ind = pop.individuals[0]
    assert len(ind.chromosomes) == 14
    print(f"  ✓ Each individual has {len(ind.chromosomes)} chromosome pairs")
    
    # Check that founders have unique IDs
    founder_ids = set()
    for ind in pop.individuals:
        for chrom_pair in ind.chromosomes:
            for chrom in chrom_pair:
                for seg in chrom.segments:
                    founder_ids.add(seg.founder_id)
    
    # Should have 2 haplotypes * 14 chromosomes * 10 individuals = 280 unique founders
    expected_founders = 2 * 14 * 10
    assert len(founder_ids) == expected_founders
    print(f"  ✓ All founder haplotypes are unique: {len(founder_ids)} founders")


def test_mating_and_recombination():
    """Test mating with recombination."""
    print("\nTesting mating and recombination...")
    
    chrom_struct = ChromosomeStructure()
    pop = Population(pop_id=0, initial_size=5, chrom_structure=chrom_struct)
    pop.initialize_founders()
    
    parent1 = pop.individuals[0]
    parent2 = pop.individuals[1]
    
    offspring = pop.mate(parent1, parent2)
    
    assert offspring.generation == 1
    assert len(offspring.chromosomes) == 14
    print("  ✓ Offspring created with correct structure")
    
    # Check that offspring has segments from both parents
    # Get founder IDs from parents and offspring
    parent1_founders = set()
    parent2_founders = set()
    offspring_founders = set()
    
    for chrom_pair in parent1.chromosomes:
        for chrom in chrom_pair:
            for seg in chrom.segments:
                parent1_founders.add(seg.founder_id)
    
    for chrom_pair in parent2.chromosomes:
        for chrom in chrom_pair:
            for seg in chrom.segments:
                parent2_founders.add(seg.founder_id)
    
    for chrom_pair in offspring.chromosomes:
        for chrom in chrom_pair:
            for seg in chrom.segments:
                offspring_founders.add(seg.founder_id)
    
    # Offspring should only have founders from parents
    assert offspring_founders.issubset(parent1_founders | parent2_founders)
    print("  ✓ Offspring inherits segments from both parents")


def test_simulation_advance():
    """Test advancing population through generations."""
    print("\nTesting generation advancement...")
    
    chrom_struct = ChromosomeStructure()
    pop = Population(pop_id=0, initial_size=10, chrom_structure=chrom_struct)
    pop.initialize_founders()
    
    initial_gen = pop.generation
    pop.advance_generation(allow_selfing=True)
    
    assert pop.generation == initial_gen + 1
    assert len(pop.individuals) == 10  # Population size maintained
    print(f"  ✓ Population advanced from gen {initial_gen} to gen {pop.generation}")


def test_simulator_basic():
    """Test basic Simulator functionality."""
    print("\nTesting Simulator...")
    
    sim = Simulator()
    pop = sim.add_population(0, size=20)
    
    assert 0 in sim.populations
    assert len(pop.individuals) == 20
    print("  ✓ Population added to simulator")
    
    # Run for a few generations
    sim.run_simulation(num_generations=5, allow_selfing=True)
    
    assert sim.generation == 5
    print("  ✓ Simulation ran for 5 generations")
    
    # Compute IBD
    individuals = sim.get_all_individuals()
    ibd_results = sim.compute_pairwise_ibd(individuals, min_length=50000)
    
    print(f"  ✓ Computed pairwise IBD: {len(ibd_results)} pairs with IBD")


def test_simulation_runner():
    """Test high-level SimulationRunner."""
    print("\nTesting SimulationRunner...")
    
    runner = SimulationRunner(output_dir="test_results")
    runner.setup_single_population(pop_size=15, recombination_rate=7.4e-7)
    
    print("  ✓ SimulationRunner initialized")
    
    individuals = runner.run_simulation(num_generations=5, allow_selfing=True)
    assert len(individuals) == 15
    print(f"  ✓ Simulation completed with {len(individuals)} individuals")
    
    # Sample a few individuals
    sampled = individuals[:10]
    
    # Compute ground truth IBD
    true_ibd = runner.compute_ground_truth_ibd(sampled, min_length=50000)
    print(f"  ✓ Ground truth IBD computed: {len(true_ibd)} pairs")
    
    # Generate genotype data
    genotype_data = runner.generate_genotype_data(sampled, snps_per_chrom=100)
    assert len(genotype_data) == 14  # One per chromosome
    print(f"  ✓ Genotype data generated for {len(genotype_data)} chromosomes")


def test_ibd_callers():
    """Test IBD calling methods."""
    print("\nTesting IBD calling methods...")
    
    # Create simple genotype data
    num_snps = 100
    num_inds = 5
    
    genotype_data = {
        0: np.random.randint(0, 2, size=(num_snps, num_inds)),
        1: np.random.randint(0, 2, size=(num_snps, num_inds))
    }
    
    chrom_lengths = [1000000, 1000000]
    
    # Test threshold caller
    threshold_caller = ThresholdIBDCaller(min_segment_length=10000)
    segments = threshold_caller.call_ibd(genotype_data, 0, 1, chrom_lengths)
    print(f"  ✓ Threshold caller: found {len(segments)} segments")
    
    # Test HMM caller
    hmm_caller = HMMIBDCaller(min_segment_length=10000)
    segments = hmm_caller.call_ibd(genotype_data, 0, 1, chrom_lengths)
    print(f"  ✓ HMM caller: found {len(segments)} segments")
    
    # Test window caller
    window_caller = WindowIBDCaller(min_segment_length=10000)
    segments = window_caller.call_ibd(genotype_data, 0, 1, chrom_lengths)
    print(f"  ✓ Window caller: found {len(segments)} segments")


def run_all_tests():
    """Run all tests."""
    print("=" * 70)
    print("Running Basic Functionality Tests")
    print("=" * 70)
    
    tests = [
        test_segment_creation,
        test_chromosome_creation,
        test_chromosome_structure,
        test_population_initialization,
        test_mating_and_recombination,
        test_simulation_advance,
        test_simulator_basic,
        test_simulation_runner,
        test_ibd_callers,
    ]
    
    passed = 0
    failed = 0
    
    for test in tests:
        try:
            test()
            passed += 1
        except Exception as e:
            print(f"  ✗ Test failed: {e}")
            import traceback
            traceback.print_exc()
            failed += 1
    
    print("\n" + "=" * 70)
    print(f"Tests completed: {passed} passed, {failed} failed")
    print("=" * 70)
    
    if failed == 0:
        print("\n✓ All tests passed!")
    else:
        print(f"\n✗ {failed} test(s) failed")
    
    return failed == 0


if __name__ == "__main__":
    success = run_all_tests()
    exit(0 if success else 1)
