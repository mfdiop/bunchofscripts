"""
IBD Calling Methods for Simulation Evaluation

Implements various IBD calling methods to compare against ground truth:
- Simple threshold-based calling
- Hidden Markov Model (HMM) approach similar to hmmIBD
- Window-based approach similar to isoRelate

Author: Generated for IBD method evaluation
Date: 2025-10-01
"""

import numpy as np
from typing import List, Tuple, Dict, Optional
from dataclasses import dataclass
from scipy.special import logsumexp
import warnings


@dataclass
class CalledIBDSegment:
    """
    Represents an IBD segment called by a detection method.
    
    Attributes:
        chrom_idx: Chromosome index
        start: Start position
        end: End position
        posterior_prob: Posterior probability of IBD (if available)
    """
    chrom_idx: int
    start: float
    end: float
    posterior_prob: float = 1.0
    
    def length(self) -> float:
        return self.end - self.start


class ThresholdIBDCaller:
    """
    Simple threshold-based IBD caller.
    
    Calls IBD when IBS (Identity by State) proportion exceeds a threshold
    in sliding windows.
    """
    
    def __init__(self, window_size: float = 50000, step_size: float = 10000,
                 ibs_threshold: float = 0.95, min_segment_length: float = 100000):
        """
        Initialize threshold caller.
        
        Args:
            window_size: Size of sliding window in bp
            step_size: Step size for sliding window in bp
            ibs_threshold: Minimum IBS proportion to call IBD
            min_segment_length: Minimum segment length to report
        """
        self.window_size = window_size
        self.step_size = step_size
        self.ibs_threshold = ibs_threshold
        self.min_segment_length = min_segment_length
    
    def call_ibd(self, genotype_data: Dict[int, np.ndarray], 
                 ind1_idx: int, ind2_idx: int,
                 chrom_lengths: List[float]) -> List[CalledIBDSegment]:
        """
        Call IBD segments between two individuals.
        
        Args:
            genotype_data: Dictionary mapping chromosome index to genotype matrix
                          (shape: num_snps x num_individuals)
            ind1_idx: Index of first individual
            ind2_idx: Index of second individual
            chrom_lengths: List of chromosome lengths
            
        Returns:
            List of called IBD segments
        """
        called_segments = []
        
        for chrom_idx, genotypes in genotype_data.items():
            if genotypes.shape[1] <= max(ind1_idx, ind2_idx):
                continue
            
            geno1 = genotypes[:, ind1_idx]
            geno2 = genotypes[:, ind2_idx]
            
            # Get SNP positions (assume evenly spaced for simplicity)
            num_snps = len(geno1)
            chrom_length = chrom_lengths[chrom_idx]
            snp_positions = np.linspace(0, chrom_length, num_snps)
            
            # Compute IBS in sliding windows
            ibd_regions = []
            current_start = None
            
            position = 0
            while position < chrom_length:
                # Find SNPs in window
                window_start = position
                window_end = min(position + self.window_size, chrom_length)
                
                in_window = (snp_positions >= window_start) & (snp_positions < window_end)
                
                if in_window.sum() > 0:
                    # Compute IBS
                    ibs = (geno1[in_window] == geno2[in_window]).mean()
                    
                    if ibs >= self.ibs_threshold:
                        if current_start is None:
                            current_start = window_start
                    else:
                        if current_start is not None:
                            # End of IBD region
                            if position - current_start >= self.min_segment_length:
                                ibd_regions.append((current_start, position))
                            current_start = None
                
                position += self.step_size
            
            # Handle final region
            if current_start is not None:
                if chrom_length - current_start >= self.min_segment_length:
                    ibd_regions.append((current_start, chrom_length))
            
            # Create called segments
            for start, end in ibd_regions:
                called_segments.append(CalledIBDSegment(
                    chrom_idx=chrom_idx,
                    start=start,
                    end=end,
                    posterior_prob=1.0
                ))
        
        return called_segments


class HMMIBDCaller:
    """
    HMM-based IBD caller similar to hmmIBD.
    
    Uses a two-state HMM (IBD vs non-IBD) with emission probabilities
    based on allele frequencies and genotype concordance.
    """
    
    def __init__(self, switch_rate: float = 1e-6, 
                 error_rate: float = 0.01,
                 min_segment_length: float = 100000,
                 min_posterior: float = 0.5):
        """
        Initialize HMM IBD caller.
        
        Args:
            switch_rate: Rate of switching between IBD/non-IBD states per bp
            error_rate: Genotyping error rate
            min_segment_length: Minimum segment length to report
            min_posterior: Minimum posterior probability to report segment
        """
        self.switch_rate = switch_rate
        self.error_rate = error_rate
        self.min_segment_length = min_segment_length
        self.min_posterior = min_posterior
    
    def compute_emission_probs(self, geno1: np.ndarray, geno2: np.ndarray,
                               allele_freqs: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """
        Compute emission probabilities for IBD and non-IBD states.
        
        Args:
            geno1: Genotypes for individual 1 (0 or 1)
            geno2: Genotypes for individual 2 (0 or 1)
            allele_freqs: Allele frequencies at each SNP
            
        Returns:
            Tuple of (p_ibd, p_not_ibd) emission probability arrays
        """
        # IBD state: genotypes should match except for errors
        match = (geno1 == geno2).astype(float)
        p_ibd = np.where(match, 1 - self.error_rate, self.error_rate)
        
        # Non-IBD state: match by chance based on allele frequencies
        # P(match | not IBD) = p^2 + (1-p)^2 for each allele
        p_match_by_chance = allele_freqs**2 + (1 - allele_freqs)**2
        p_not_ibd = np.where(match, p_match_by_chance, 1 - p_match_by_chance)
        
        # Avoid log(0) issues
        p_ibd = np.clip(p_ibd, 1e-10, 1 - 1e-10)
        p_not_ibd = np.clip(p_not_ibd, 1e-10, 1 - 1e-10)
        
        return p_ibd, p_not_ibd
    
    def forward_backward(self, emission_ibd: np.ndarray, emission_not_ibd: np.ndarray,
                        positions: np.ndarray) -> np.ndarray:
        """
        Run forward-backward algorithm to compute posterior probabilities.
        
        Args:
            emission_ibd: Emission probabilities for IBD state
            emission_not_ibd: Emission probabilities for non-IBD state
            positions: Physical positions of SNPs
            
        Returns:
            Posterior probability of IBD at each SNP
        """
        n = len(emission_ibd)
        
        # Log space to avoid underflow
        log_emission_ibd = np.log(emission_ibd)
        log_emission_not_ibd = np.log(emission_not_ibd)
        
        # Compute transition probabilities based on distance
        distances = np.diff(positions)
        distances = np.clip(distances, 0, None)
        
        # P(switch) = 1 - exp(-switch_rate * distance)
        p_switch = 1 - np.exp(-self.switch_rate * distances)
        p_stay = 1 - p_switch
        
        # Initialize forward variables (log space)
        # Start with equal prior
        log_alpha = np.zeros((n, 2))  # [IBD, not IBD]
        log_alpha[0] = [np.log(0.5) + log_emission_ibd[0],
                       np.log(0.5) + log_emission_not_ibd[0]]
        
        # Forward pass
        for t in range(1, n):
            # Transition matrix (log space)
            log_trans = np.log([[p_stay[t-1], p_switch[t-1]],
                               [p_switch[t-1], p_stay[t-1]]])
            
            # Update
            for s in range(2):  # Current state
                log_alpha[t, s] = logsumexp(log_alpha[t-1] + log_trans[:, s])
            
            # Add emissions
            log_alpha[t, 0] += log_emission_ibd[t]
            log_alpha[t, 1] += log_emission_not_ibd[t]
        
        # Initialize backward variables
        log_beta = np.zeros((n, 2))
        log_beta[-1] = [0, 0]  # log(1) = 0
        
        # Backward pass
        for t in range(n - 2, -1, -1):
            log_trans = np.log([[p_stay[t], p_switch[t]],
                               [p_switch[t], p_stay[t]]])
            
            for s in range(2):  # Current state
                emissions = np.array([log_emission_ibd[t+1], log_emission_not_ibd[t+1]])
                log_beta[t, s] = logsumexp(log_trans[s, :] + emissions + log_beta[t+1])
        
        # Compute posteriors
        log_gamma = log_alpha + log_beta
        log_gamma -= logsumexp(log_gamma, axis=1, keepdims=True)  # Normalize
        
        # Return posterior probability of IBD state
        return np.exp(log_gamma[:, 0])
    
    def call_ibd(self, genotype_data: Dict[int, np.ndarray],
                 ind1_idx: int, ind2_idx: int,
                 chrom_lengths: List[float],
                 allele_freqs: Optional[Dict[int, np.ndarray]] = None) -> List[CalledIBDSegment]:
        """
        Call IBD segments using HMM approach.
        
        Args:
            genotype_data: Dictionary mapping chromosome index to genotype matrix
            ind1_idx: Index of first individual
            ind2_idx: Index of second individual
            chrom_lengths: List of chromosome lengths
            allele_freqs: Optional allele frequencies (computed from data if not provided)
            
        Returns:
            List of called IBD segments
        """
        called_segments = []
        
        for chrom_idx, genotypes in genotype_data.items():
            if genotypes.shape[1] <= max(ind1_idx, ind2_idx):
                continue
            
            geno1 = genotypes[:, ind1_idx]
            geno2 = genotypes[:, ind2_idx]
            
            # Compute allele frequencies if not provided
            if allele_freqs is None or chrom_idx not in allele_freqs:
                afreqs = genotypes.mean(axis=1)
                afreqs = np.clip(afreqs, 0.01, 0.99)  # Avoid extreme frequencies
            else:
                afreqs = allele_freqs[chrom_idx]
            
            # Get SNP positions
            num_snps = len(geno1)
            chrom_length = chrom_lengths[chrom_idx]
            snp_positions = np.linspace(0, chrom_length, num_snps)
            
            # Compute emission probabilities
            p_ibd, p_not_ibd = self.compute_emission_probs(geno1, geno2, afreqs)
            
            # Run forward-backward
            posterior_ibd = self.forward_backward(p_ibd, p_not_ibd, snp_positions)
            
            # Call segments where posterior exceeds threshold
            in_segment = posterior_ibd >= self.min_posterior
            
            # Find contiguous regions
            segment_starts = []
            segment_ends = []
            
            if in_segment[0]:
                segment_starts.append(0)
            
            for i in range(1, len(in_segment)):
                if in_segment[i] and not in_segment[i-1]:
                    segment_starts.append(i)
                elif not in_segment[i] and in_segment[i-1]:
                    segment_ends.append(i)
            
            if in_segment[-1]:
                segment_ends.append(len(in_segment))
            
            # Create segments
            for start_idx, end_idx in zip(segment_starts, segment_ends):
                start_pos = snp_positions[start_idx]
                end_pos = snp_positions[end_idx - 1] if end_idx > 0 else chrom_length
                
                if end_pos - start_pos >= self.min_segment_length:
                    # Compute mean posterior probability
                    mean_posterior = posterior_ibd[start_idx:end_idx].mean()
                    
                    called_segments.append(CalledIBDSegment(
                        chrom_idx=chrom_idx,
                        start=start_pos,
                        end=end_pos,
                        posterior_prob=mean_posterior
                    ))
        
        return called_segments


class WindowIBDCaller:
    """
    Window-based IBD caller similar to isoRelate approach.
    
    Uses overlapping windows and requires a minimum number of informative sites.
    """
    
    def __init__(self, window_size: float = 100000,
                 overlap: float = 50000,
                 min_informative_sites: int = 5,
                 ibs_threshold: float = 0.9,
                 min_segment_length: float = 100000):
        """
        Initialize window-based caller.
        
        Args:
            window_size: Window size in bp
            overlap: Overlap between windows in bp
            min_informative_sites: Minimum number of polymorphic sites in window
            ibs_threshold: IBS threshold for calling IBD
            min_segment_length: Minimum segment length to report
        """
        self.window_size = window_size
        self.overlap = overlap
        self.min_informative_sites = min_informative_sites
        self.ibs_threshold = ibs_threshold
        self.min_segment_length = min_segment_length
    
    def call_ibd(self, genotype_data: Dict[int, np.ndarray],
                 ind1_idx: int, ind2_idx: int,
                 chrom_lengths: List[float]) -> List[CalledIBDSegment]:
        """
        Call IBD segments using window-based approach.
        
        Args:
            genotype_data: Dictionary mapping chromosome index to genotype matrix
            ind1_idx: Index of first individual
            ind2_idx: Index of second individual
            chrom_lengths: List of chromosome lengths
            
        Returns:
            List of called IBD segments
        """
        called_segments = []
        
        for chrom_idx, genotypes in genotype_data.items():
            if genotypes.shape[1] <= max(ind1_idx, ind2_idx):
                continue
            
            geno1 = genotypes[:, ind1_idx]
            geno2 = genotypes[:, ind2_idx]
            
            # Identify polymorphic sites
            allele_freqs = genotypes.mean(axis=1)
            is_polymorphic = (allele_freqs > 0.01) & (allele_freqs < 0.99)
            
            # Get SNP positions
            num_snps = len(geno1)
            chrom_length = chrom_lengths[chrom_idx]
            snp_positions = np.linspace(0, chrom_length, num_snps)
            
            # Slide windows
            step = self.window_size - self.overlap
            position = 0
            window_ibd_status = []
            window_positions = []
            
            while position < chrom_length:
                window_end = min(position + self.window_size, chrom_length)
                
                # Find SNPs in window
                in_window = (snp_positions >= position) & (snp_positions < window_end)
                in_window_polymorphic = in_window & is_polymorphic
                
                num_informative = in_window_polymorphic.sum()
                
                if num_informative >= self.min_informative_sites:
                    # Compute IBS at informative sites
                    ibs = (geno1[in_window_polymorphic] == 
                          geno2[in_window_polymorphic]).mean()
                    
                    is_ibd = ibs >= self.ibs_threshold
                else:
                    is_ibd = False
                
                window_ibd_status.append(is_ibd)
                window_positions.append((position, window_end))
                
                position += step
            
            # Merge adjacent IBD windows
            if not window_ibd_status:
                continue
            
            segment_start = None
            for i, is_ibd in enumerate(window_ibd_status):
                if is_ibd:
                    if segment_start is None:
                        segment_start = window_positions[i][0]
                else:
                    if segment_start is not None:
                        segment_end = window_positions[i-1][1]
                        if segment_end - segment_start >= self.min_segment_length:
                            called_segments.append(CalledIBDSegment(
                                chrom_idx=chrom_idx,
                                start=segment_start,
                                end=segment_end,
                                posterior_prob=1.0
                            ))
                        segment_start = None
            
            # Handle final segment
            if segment_start is not None:
                segment_end = window_positions[-1][1]
                if segment_end - segment_start >= self.min_segment_length:
                    called_segments.append(CalledIBDSegment(
                        chrom_idx=chrom_idx,
                        start=segment_start,
                        end=segment_end,
                        posterior_prob=1.0
                    ))
        
        return called_segments


def evaluate_ibd_calling(true_ibd: List[Tuple[int, float, float]],
                        called_ibd: List[CalledIBDSegment],
                        total_genome_length: float) -> Dict[str, float]:
    """
    Evaluate IBD calling performance against ground truth.
    
    Args:
        true_ibd: Ground truth IBD segments as (chrom_idx, start, end)
        called_ibd: Called IBD segments
        total_genome_length: Total genome length for computing rates
        
    Returns:
        Dictionary with performance metrics (sensitivity, precision, etc.)
    """
    # Convert to sets of covered base pairs
    true_covered = set()
    for chrom_idx, start, end in true_ibd:
        for pos in range(int(start), int(end)):
            true_covered.add((chrom_idx, pos))
    
    called_covered = set()
    for seg in called_ibd:
        for pos in range(int(seg.start), int(seg.end)):
            called_covered.add((seg.chrom_idx, pos))
    
    # Compute metrics
    true_positive = len(true_covered & called_covered)
    false_positive = len(called_covered - true_covered)
    false_negative = len(true_covered - called_covered)
    
    sensitivity = true_positive / len(true_covered) if len(true_covered) > 0 else 0
    precision = true_positive / len(called_covered) if len(called_covered) > 0 else 0
    f1_score = (2 * sensitivity * precision / (sensitivity + precision) 
                if sensitivity + precision > 0 else 0)
    
    return {
        'sensitivity': sensitivity,
        'precision': precision,
        'f1_score': f1_score,
        'true_positive_bp': true_positive,
        'false_positive_bp': false_positive,
        'false_negative_bp': false_negative,
    }
