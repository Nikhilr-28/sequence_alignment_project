import sys
import time
import tracemalloc

def generate_string(base_string, indices):
    """
    Generate string according to the insertion rules
    At each step, insert a copy of the current string after the specified position
    """
    current = base_string
    
    for index in indices:
        # Insert a copy of current string after position 'index'
        if index >= len(current):
            # If index is beyond string length, append at end
            current = current + current
        else:
            # Insert copy of current string after position 'index'
            current = current[:index+1] + current + current[index+1:]
    
    return current

def read_input_file(filename):
    """
    Read input file and return parsed data
    Better parsing that handles identical s1 and s2
    """
    with open(filename, 'r') as f:
        lines = [line.strip() for line in f.readlines()]
    
    # Find all DNA strings (lines that contain only A, C, G, T and are reasonable length)
    dna_lines = []
    for i, line in enumerate(lines):
        if all(c in 'ACGT' for c in line) and len(line) >= 3:
            dna_lines.append((i, line))
    
    # If we found exactly 2 DNA strings, use those as s1 and s2
    if len(dna_lines) == 2:
        s1_line, s1 = dna_lines[0]
        s2_line, s2 = dna_lines[1]
        
        # All numbers between s1 and s2 are indices for s1
        s1_indices = []
        for i in range(s1_line + 1, s2_line):
            try:
                idx = int(lines[i])
                s1_indices.append(idx)
            except:
                pass
        
        # All numbers after s2 are indices for s2
        s2_indices = []
        for i in range(s2_line + 1, len(lines)):
            try:
                idx = int(lines[i])
                s2_indices.append(idx)
            except:
                pass
    
    # Fallback to old logic if needed
    else:
        s1 = lines[0]
        s1_indices = []
        s2_start = -1
        
        # Find where s2 starts (look for another DNA string)
        for i in range(1, len(lines)):
            line = lines[i]
            # Check if it's a DNA string
            if all(c in 'ACGT' for c in line) and len(line) >= 3 and line != s1:
                s2_start = i
                break
            else:
                # Try to convert to int
                try:
                    idx = int(line)
                    s1_indices.append(idx)
                except:
                    pass
        
        # Get s2 and its indices
        s2 = lines[s2_start]
        s2_indices = []
        for i in range(s2_start + 1, len(lines)):
            try:
                idx = int(lines[i])
                s2_indices.append(idx)
            except:
                pass
    
    return s1, s1_indices, s2, s2_indices

# Mismatch costs table
alpha = {
    ('A', 'A'): 0, ('A', 'C'): 110, ('A', 'G'): 48, ('A', 'T'): 94,
    ('C', 'A'): 110, ('C', 'C'): 0, ('C', 'G'): 118, ('C', 'T'): 48,
    ('G', 'A'): 48, ('G', 'C'): 118, ('G', 'G'): 0, ('G', 'T'): 110,
    ('T', 'A'): 94, ('T', 'C'): 48, ('T', 'G'): 110, ('T', 'T'): 0
}

# Gap penalty
delta = 30

def sequence_alignment_efficient(s1, s2):
    """
    Memory-efficient sequence alignment using divide and conquer
    Uses O(min(m,n)) space instead of O(m*n)
    """
    m = len(s1)
    n = len(s2)
    
    # If one string is empty, just add gaps
    if m == 0:
        return n * delta, "_" * n, s2
    if n == 0:
        return m * delta, s1, "_" * m
    
    # If one string is very short, use basic algorithm
    if min(m, n) <= 3:
        return sequence_alignment_basic(s1, s2)
    
    # Divide and conquer
    mid = m // 2
    
    # Compute optimal cost for left half
    left_costs = compute_costs_forward(s1[:mid], s2)
    
    # Compute optimal cost for right half (in reverse)
    right_costs = compute_costs_backward(s1[mid:], s2)
    
    # Find optimal split point
    total_costs = [left + right for left, right in zip(left_costs, right_costs)]
    min_cost = min(total_costs)
    split_j = total_costs.index(min_cost)
    
    # Recursively solve subproblems
    left_cost, left_align1, left_align2 = sequence_alignment_efficient(s1[:mid], s2[:split_j])
    right_cost, right_align1, right_align2 = sequence_alignment_efficient(s1[mid:], s2[split_j:])
    
    # Combine results
    total_cost = left_cost + right_cost
    aligned_s1 = left_align1 + right_align1
    aligned_s2 = left_align2 + right_align2
    
    return total_cost, aligned_s1, aligned_s2

def compute_costs_forward(s1, s2):
    """
    Compute alignment costs from left to right, keeping only current and previous columns
    """
    m = len(s1)
    n = len(s2)
    
    # Only keep two columns
    prev_col = [j * delta for j in range(n + 1)]
    curr_col = [0] * (n + 1)
    
    for i in range(1, m + 1):
        curr_col[0] = i * delta
        for j in range(1, n + 1):
            match_cost = prev_col[j-1] + alpha[(s1[i-1], s2[j-1])]
            gap_s2 = prev_col[j] + delta
            gap_s1 = curr_col[j-1] + delta
            curr_col[j] = min(match_cost, gap_s2, gap_s1)
        prev_col, curr_col = curr_col, prev_col
    
    return prev_col

def compute_costs_backward(s1, s2):
    """
    Compute alignment costs from right to left, keeping only current and previous columns
    """
    m = len(s1)
    n = len(s2)
    
    # Only keep two columns
    prev_col = [j * delta for j in range(n, -1, -1)]
    curr_col = [0] * (n + 1)
    
    for i in range(m - 1, -1, -1):
        curr_col[n] = (m - i) * delta
        for j in range(n - 1, -1, -1):
            match_cost = prev_col[j+1] + alpha[(s1[i], s2[j])]
            gap_s2 = prev_col[j] + delta
            gap_s1 = curr_col[j+1] + delta
            curr_col[j] = min(match_cost, gap_s2, gap_s1)
        prev_col, curr_col = curr_col, prev_col
    
    return prev_col

def sequence_alignment_basic(s1, s2):
    """
    Basic DP algorithm for small cases
    """
    m = len(s1)
    n = len(s2)
    
    # Create DP table
    dp = [[0 for _ in range(n + 1)] for _ in range(m + 1)]
    
    # Initialize base cases
    for i in range(m + 1):
        dp[i][0] = i * delta
    for j in range(n + 1):
        dp[0][j] = j * delta
    
    # Fill the DP table
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            match_cost = dp[i-1][j-1] + alpha[(s1[i-1], s2[j-1])]
            gap_s2 = dp[i-1][j] + delta
            gap_s1 = dp[i][j-1] + delta
            dp[i][j] = min(match_cost, gap_s2, gap_s1)
    
    # Traceback
    aligned_s1 = ""
    aligned_s2 = ""
    i, j = m, n
    
    while i > 0 or j > 0:
        if i > 0 and j > 0:
            if dp[i][j] == dp[i-1][j-1] + alpha[(s1[i-1], s2[j-1])]:
                aligned_s1 = s1[i-1] + aligned_s1
                aligned_s2 = s2[j-1] + aligned_s2
                i -= 1
                j -= 1
            elif dp[i][j] == dp[i-1][j] + delta:
                aligned_s1 = s1[i-1] + aligned_s1
                aligned_s2 = "_" + aligned_s2
                i -= 1
            else:
                aligned_s1 = "_" + aligned_s1
                aligned_s2 = s2[j-1] + aligned_s2
                j -= 1
        elif i > 0:
            aligned_s1 = s1[i-1] + aligned_s1
            aligned_s2 = "_" + aligned_s2
            i -= 1
        else:
            aligned_s1 = "_" + aligned_s1
            aligned_s2 = s2[j-1] + aligned_s2
            j -= 1
    
    return dp[m][n], aligned_s1, aligned_s2

def measure_time_and_memory(func, *args):
    """
    Measure execution time and memory usage of a function
    """
    # Start memory tracking
    tracemalloc.start()
    start_time = time.time()
    
    # Execute function
    result = func(*args)
    
    # Measure time and memory
    end_time = time.time()
    current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    
    execution_time = (end_time - start_time) * 1000  # Convert to milliseconds
    memory_usage = peak / 1024  # Convert to KB
    
    return result, execution_time, memory_usage

def main():
    if len(sys.argv) != 3:
        print("Usage: python efficient_3.py input_file output_file")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    # Read input file
    s1, s1_indices, s2, s2_indices = read_input_file(input_file)
    
    # Generate strings
    gen_s1 = generate_string(s1, s1_indices)
    gen_s2 = generate_string(s2, s2_indices)
    
    # Perform sequence alignment with time and memory measurement
    def alignment_func():
        return sequence_alignment_efficient(gen_s1, gen_s2)
    
    result, exec_time, memory_usage = measure_time_and_memory(alignment_func)
    cost, aligned_s1, aligned_s2 = result
    
    # Write output
    with open(output_file, 'w') as f:
        f.write(f"{cost}\n")
        f.write(f"{aligned_s1}\n")
        f.write(f"{aligned_s2}\n")
        f.write(f"{exec_time}\n")
        f.write(f"{memory_usage}\n")

if __name__ == "__main__":
    main()