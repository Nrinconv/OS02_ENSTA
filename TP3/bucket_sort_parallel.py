from mpi4py import MPI
import numpy as np
import time

def bucket_sort(arr):
    #bucket sort implementation
    if len(arr) == 0:
        return arr  

    min_val, max_val = min(arr), max(arr)
    bucket_count = len(arr)
    # init buckets
    buckets = [[] for _ in range(bucket_count)]

    for num in arr:
        index = int(bucket_count * (num - min_val) / (max_val - min_val + 1e-6))
        buckets[index].append(num)

    # sort each bucket
    sorted_array = []
    for bucket in buckets:
        sorted_array.extend(sorted(bucket))

    return sorted_array

# init MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()  
size = comm.Get_size()

# numbers of elements to sort
N = 1000000

start_time = time.time()

if rank == 0:
    # random numbers
    data = np.random.rand(N)  
    #print(f"Original data: {data}")

    # find optimal bucket intervals based on actual data distribution (the reason for the percentile)
    bucket_intervals = np.percentile(data, np.linspace(0, 100, size + 1))

    # distribute numbers to the correct bucket
    local_data = [[] for _ in range(size)]
    for num in data:
        for i in range(size):
            if bucket_intervals[i] <= num < bucket_intervals[i + 1]:
                local_data[i].append(num)
                break

    for i in range(1, size):
        comm.send(local_data[i], dest=i, tag=0)

    local_bucket = local_data[0]

else:
    local_bucket = comm.recv(source=0, tag=0)

# each process sorts its own bucket
sorted_local_bucket = bucket_sort(local_bucket)

# gather sorted buckets in process 0
sorted_data = comm.gather(sorted_local_bucket, root=0)

# process 0 merges all sorted lists
if rank == 0:
    final_sorted_data = np.concatenate(sorted_data)  # Flatten the list
    #print(f"\nFinal sorted data: {final_sorted_data}")

    # Stop measuring total execution time and print the result
    end_time = time.time()
    print(f"\nExecution time: {end_time - start_time:.6f}")
