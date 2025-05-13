import subprocess
import os
import pandas as pd
import matplotlib.pyplot as plt

# Create results directories
os.makedirs('results_basic', exist_ok=True)
os.makedirs('results_efficient', exist_ok=True)

print("Running all 15 datapoint tests...")
results = []

# Run tests on all 15 datapoint files
for i in range(1, 16):
    input_file = f'datapoints/in{i}.txt'
    basic_output = f'results_basic/result{i}.txt'
    efficient_output = f'results_efficient/result{i}.txt'
    
    print(f"Processing in{i}.txt...")
    
    # Run basic algorithm
    subprocess.run(['python', 'basic_3.py', input_file, basic_output])
    
    # Run efficient algorithm
    subprocess.run(['python', 'efficient_3.py', input_file, efficient_output])
    
    # Read results and calculate problem size
    try:
        # Read basic results
        with open(basic_output, 'r') as f:
            basic_lines = f.readlines()
            basic_cost = basic_lines[0].strip()
            basic_time = float(basic_lines[3].strip())
            basic_memory = float(basic_lines[4].strip())
        
        # Read efficient results  
        with open(efficient_output, 'r') as f:
            efficient_lines = f.readlines()
            efficient_cost = efficient_lines[0].strip()
            efficient_time = float(efficient_lines[3].strip())
            efficient_memory = float(efficient_lines[4].strip())
        
        # Calculate problem size (M+N)
        from basic_3 import read_input_file, generate_string
        s1, s1_indices, s2, s2_indices = read_input_file(input_file)
        gen_s1 = generate_string(s1, s1_indices)
        gen_s2 = generate_string(s2, s2_indices)
        problem_size = len(gen_s1) + len(gen_s2)
        
        result = {
            'M+N': problem_size,
            'Time (Basic)': round(basic_time, 3),
            'Time (Efficient)': round(efficient_time, 3),
            'Memory (Basic)': round(basic_memory, 3),
            'Memory (Efficient)': round(efficient_memory, 3)
        }
        results.append(result)
        
        print(f"  M+N={problem_size}, Basic: {basic_time:.3f}ms/{basic_memory:.3f}KB, Efficient: {efficient_time:.3f}ms/{efficient_memory:.3f}KB")
    
    except Exception as e:
        print(f"Error processing in{i}.txt: {e}")

# Sort results by problem size
results = sorted(results, key=lambda x: x['M+N'])

# Create DataFrame for easier handling
df = pd.DataFrame(results)

# Create the plots
plt.style.use('default')
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 7))

# Graph 1: Memory vs Problem Size
ax1.plot(df['M+N'], df['Memory (Basic)'], 'bo-', label='Basic Algorithm', linewidth=3, markersize=8)
ax1.plot(df['M+N'], df['Memory (Efficient)'], 'ro-', label='Memory Efficient Algorithm', linewidth=3, markersize=8)
ax1.set_xlabel('Problem Size (M+N)', fontsize=14, fontweight='bold')
ax1.set_ylabel('Memory (KB)', fontsize=14, fontweight='bold')
ax1.set_title('Memory Usage vs Problem Size', fontsize=16, fontweight='bold')
ax1.legend(fontsize=12)
ax1.grid(True, alpha=0.3)
ax1.set_xlim(0, max(df['M+N']) * 1.05)

# Graph 2: Time vs Problem Size  
ax2.plot(df['M+N'], df['Time (Basic)'], 'bo-', label='Basic Algorithm', linewidth=3, markersize=8)
ax2.plot(df['M+N'], df['Time (Efficient)'], 'ro-', label='Memory Efficient Algorithm', linewidth=3, markersize=8)
ax2.set_xlabel('Problem Size (M+N)', fontsize=14, fontweight='bold')
ax2.set_ylabel('Time (milliseconds)', fontsize=14, fontweight='bold')
ax2.set_title('CPU Time vs Problem Size', fontsize=16, fontweight='bold')
ax2.legend(fontsize=12)
ax2.grid(True, alpha=0.3)
ax2.set_xlim(0, max(df['M+N']) * 1.05)

# Save individual graphs for Summary.docx
plt.tight_layout()
plt.savefig('combined_plots.png', dpi=300, bbox_inches='tight')

# Save individual graphs separately
fig1, ax1_alone = plt.subplots(figsize=(10, 7))
ax1_alone.plot(df['M+N'], df['Memory (Basic)'], 'bo-', label='Basic Algorithm', linewidth=3, markersize=8)
ax1_alone.plot(df['M+N'], df['Memory (Efficient)'], 'ro-', label='Memory Efficient Algorithm', linewidth=3, markersize=8)
ax1_alone.set_xlabel('Problem Size (M+N)', fontsize=14, fontweight='bold')
ax1_alone.set_ylabel('Memory (KB)', fontsize=14, fontweight='bold')
ax1_alone.set_title('Memory Usage vs Problem Size', fontsize=16, fontweight='bold')
ax1_alone.legend(fontsize=12)
ax1_alone.grid(True, alpha=0.3)
ax1_alone.set_xlim(0, max(df['M+N']) * 1.05)
plt.tight_layout()
plt.savefig('Graph1_Memory_vs_ProblemSize.png', dpi=300, bbox_inches='tight')
plt.close()

fig2, ax2_alone = plt.subplots(figsize=(10, 7))
ax2_alone.plot(df['M+N'], df['Time (Basic)'], 'bo-', label='Basic Algorithm', linewidth=3, markersize=8)
ax2_alone.plot(df['M+N'], df['Time (Efficient)'], 'ro-', label='Memory Efficient Algorithm', linewidth=3, markersize=8)
ax2_alone.set_xlabel('Problem Size (M+N)', fontsize=14, fontweight='bold')
ax2_alone.set_ylabel('Time (milliseconds)', fontsize=14, fontweight='bold')
ax2_alone.set_title('CPU Time vs Problem Size', fontsize=16, fontweight='bold')
ax2_alone.legend(fontsize=12)
ax2_alone.grid(True, alpha=0.3)
ax2_alone.set_xlim(0, max(df['M+N']) * 1.05)
plt.tight_layout()
plt.savefig('Graph2_Time_vs_ProblemSize.png', dpi=300, bbox_inches='tight')
plt.close()

# Print table data for Summary.docx
print("\n=== TABLE DATA FOR SUMMARY.DOCX ===")
print("M+N\t\tTime (Basic)\tTime (Efficient)\tMemory (Basic)\tMemory (Efficient)")
print("-" * 80)
for _, row in df.iterrows():
    print(f"{row['M+N']}\t\t{row['Time (Basic)']}\t\t{row['Time (Efficient)']}\t\t{row['Memory (Basic)']}\t\t{row['Memory (Efficient)']}")

# Analysis insights
print("\n=== INSIGHTS FOR SUMMARY.DOCX ===")
print("\nGraph1 -- Memory vs Problem Size (M+N)")
print("Nature of the Graph:")
print("Basic: Polynomial (Quadratic)")
print("Efficient: Linear")
print("\nExplanation:")
print("The basic algorithm shows quadratic memory growth due to its O(m×n) space complexity, ")
print("while the memory-efficient algorithm maintains linear growth with O(min(m,n)) space usage.")

print("\nGraph2 -- Time vs Problem Size (M+N)")
print("Nature of the Graph:")
print("Basic: Polynomial (Quadratic)")  
print("Efficient: Polynomial (Quadratic)")
print("\nExplanation:")
print("Both algorithms exhibit quadratic time complexity O(m×n), though the memory-efficient ")
print("algorithm shows slightly higher execution times due to recursive overhead from the divide-and-conquer approach.")

print("\n=== FILES CREATED ===")
print("1. Graph1_Memory_vs_ProblemSize.png - Insert this in Graph1 section")
print("2. Graph2_Time_vs_ProblemSize.png - Insert this in Graph2 section")
print("3. combined_plots.png - Shows both graphs together")
print("4. Table data printed above - Copy to Summary.docx")