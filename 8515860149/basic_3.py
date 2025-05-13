import sys
import time
import tracemalloc

def generate_string(base_string, indices):
    """Generate string by inserting copies at specified positions"""
    current = base_string
    
    for index in indices:
        if index >= len(current):
            current = current + current
        else:
            current = current[:index+1] + current + current[index+1:]
    
    return current

def read_input_file(filename):
    """Parse input file to extract base strings and indices"""
    with open(filename, 'r') as f:
        lines = [line.strip() for line in f.readlines()]
    
    # Find DNA strings in the file
    dna_lines = []
    for i, line in enumerate(lines):
        if all(c in 'ACGT' for c in line) and len(line) >= 3:
            dna_lines.append((i, line))
    
    # Extract indices based on DNA string positions
    if len(dna_lines) == 2:
        s1_line, s1 = dna_lines[0]
        s2_line, s2 = dna_lines[1]
        
        s1_indices = []
        for i in range(s1_line + 1, s2_line):
            try:
                s1_indices.append(int(lines[i]))
            except:
                pass
        
        s2_indices = []
        for i in range(s2_line + 1, len(lines)):
            try:
                s2_indices.append(int(lines[i]))
            except:
                pass
    else:
        # Fallback parsing method
        s1 = lines[0]
        s1_indices = []
        s2_start = -1
        
        for i in range(1, len(lines)):
            line = lines[i]
            if all(c in 'ACGT' for c in line) and len(line) >= 3 and line != s1:
                s2_start = i
                break
            else:
                try:
                    s1_indices.append(int(line))
                except:
                    pass
        
        s2 = lines[s2_start]
        s2_indices = []
        for i in range(s2_start + 1, len(lines)):
            try:
                s2_indices.append(int(lines[i]))
            except:
                pass
    
    return s1, s1_indices, s2, s2_indices

# Mismatch cost matrix
alpha = {
    ('A', 'A'): 0, ('A', 'C'): 110, ('A', 'G'): 48, ('A', 'T'): 94,
    ('C', 'A'): 110, ('C', 'C'): 0, ('C', 'G'): 118, ('C', 'T'): 48,
    ('G', 'A'): 48, ('G', 'C'): 118, ('G', 'G'): 0, ('G', 'T'): 110,
    ('T', 'A'): 94, ('T', 'C'): 48, ('T', 'G'): 110, ('T', 'T'): 0
}

delta = 30  # Gap penalty

def sequence_alignment(s1, s2):
    """Basic DP algorithm for sequence alignment"""
    m, n = len(s1), len(s2)
    
    # Initialize DP table
    dp = [[0 for _ in range(n + 1)] for _ in range(m + 1)]
    
    for i in range(m + 1):
        dp[i][0] = i * delta
    for j in range(n + 1):
        dp[0][j] = j * delta
    
    # Fill DP table
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            match_cost = dp[i-1][j-1] + alpha[(s1[i-1], s2[j-1])]
            gap_s2 = dp[i-1][j] + delta
            gap_s1 = dp[i][j-1] + delta
            dp[i][j] = min(match_cost, gap_s2, gap_s1)
    
    # Traceback to reconstruct alignment
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
    """Measure function execution time and memory usage"""
    tracemalloc.start()
    start_time = time.time()
    
    result = func(*args)
    
    end_time = time.time()
    current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    
    execution_time = (end_time - start_time) * 1000
    memory_usage = peak / 1024
    
    return result, execution_time, memory_usage

def main():
    if len(sys.argv) != 3:
        print("Usage: python basic_3.py input_file output_file")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    # Parse input and generate strings
    s1, s1_indices, s2, s2_indices = read_input_file(input_file)
    gen_s1 = generate_string(s1, s1_indices)
    gen_s2 = generate_string(s2, s2_indices)
    
    # Run alignment with performance measurement
    def alignment_func():
        return sequence_alignment(gen_s1, gen_s2)
    
    result, exec_time, memory_usage = measure_time_and_memory(alignment_func)
    cost, aligned_s1, aligned_s2 = result
    
    # Write results to output file
    with open(output_file, 'w') as f:
        f.write(f"{cost}\n")
        f.write(f"{aligned_s1}\n")
        f.write(f"{aligned_s2}\n")
        f.write(f"{exec_time}\n")
        f.write(f"{memory_usage}\n")

if __name__ == "__main__":
    main()