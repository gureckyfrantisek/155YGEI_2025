import random
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
import numpy as np
# K-Means clustering

def gen_pts(dimensions, cluster_count):
    points = []
    
    for i in range(cluster_count):
        # Random sized clusters
        cluster_points = random.randint(10, 200)
        
        # Random offset from center
        offset = [0] * dimensions
        for o in range(dimensions):
            offset[o] = (random.random() - 0.5) * 50        # Random float from -25 to 25
        
        # Generate cluster points
        for j in range(cluster_points):
            point = []
            for k in range(dimensions):
                point.append((random.random() - 0.5) * 10 + offset[k])   # Random float from -5 to 5 with offset
            
            points.append(point)
                
    return points

def distance(point_a, point_b):
    dim = len(point_a)
    sum_delta_sq = 0
    
    # Calculate the sum of deltas squared
    for i in range(dim):
        sum_delta_sq += (point_a[i] - point_b[i]) ** 2
    
    # The distance is the square root
    return sum_delta_sq ** 0.5

def k_means(points, cluster_count):
    print("Running k-means")
    dim = len(points[0])
    point_count = len(points)
    
    centroids = []
    
    # Generate centroids randomly
    for i in range(cluster_count):
        centroid = []
        for j in range(dim):
            centroid.append((random.random() - 0.5) * 40)   # Random float from -20 to 20 to match the points generation
        
        centroids.append(centroid)
    
    print(centroids)
    # Initialize the iteration
    done = False
    belongs_to = [None] * point_count
    iteration_count = 0
    
    while not done:
        iteration_count += 1
        print(f"Starting {iteration_count}. iteration")
        done = True
        
        # Storing the old belongs_to for comparison as the terminal condition
        old_belongs_to = belongs_to.copy()
        
        # Find closest cluster for all points
        for i in range(point_count):
            nearest = 9999999
            for j, center in enumerate(centroids):
                dist = distance(points[i], center)
                if dist < nearest:
                    belongs_to[i] = j
                    nearest = dist
        
        # Recalculate the centroids
        for i, center in enumerate(centroids):
            new_center = [0] * dim
            children = 0
            
            for j in range(point_count):
                # If this point belongs to this centroid, add it to the sum
                if i == belongs_to[j]:
                    children += 1
                    
                    for k in range(dim):
                        new_center[k] += points[j][k]
            
            # We've passed all the points for the centroid, divide the coordinates by the belonging count
            for k in range(dim):
                new_center[k] /= max(children, 1)
            
            # Update the centroids
            centroids[i] = new_center
        
        if old_belongs_to != belongs_to:
            done = False
    
    # After the iteration is over, return the belonging array and the final centroids
    return belongs_to, centroids

def visualize(points, belongs_to, centroids):
    # Get the same coordinates to list, if the dimension is two
    if len(points[0]) > 2:
        print("Unable to visualize in more than 2D")
        return
    
    cluster_count = len(centroids)
    
    # Plot points for the cluster
    for i in range(cluster_count):
        # Separate lists for x and y
        cluster_x = []
        cluster_y = []
        
        for j in range(len(points)):
            if belongs_to[j] == i:
                cluster_x.append(points[j][0])
                cluster_y.append(points[j][1])

        plt.scatter(cluster_x, cluster_y, s=10)

    # Plot centroids
    cent_x = [c[0] for c in centroids]
    cent_y = [c[1] for c in centroids]
    plt.scatter(cent_x, cent_y, color="black", s=200, marker="X")

    plt.title("K-Means Clustering Visualization")
    plt.show()
    
def k_means_scikit(cluster_count, init, points):
    # First translate the point data to np.array
    data = np.array(points)
    kmeans_instance = KMeans(n_clusters=cluster_count, random_state=0, n_init="auto", init=init)
    kmeans = kmeans_instance.fit(data)
    
    # Prepare the data as simple lists and return them
    belongs_to = kmeans.labels_.tolist()               # Similar to our belongs_to array of indices
    centroids = kmeans.cluster_centers_.tolist()     # The centroids of clusters

    return belongs_to, centroids
    
def main():
    dim = 2
    cluster_count = 5
    
    points = gen_pts(dim, cluster_count)
    
    belongs_to, centroids = k_means(points, cluster_count)
    belongs_to_sc, centroids_sc = k_means_scikit(cluster_count, "random", points)
    
    # Vizualize the results
    visualize(points, belongs_to, centroids)
    visualize(points, belongs_to_sc, centroids_sc)

if __name__ == "__main__":
    main()