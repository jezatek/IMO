import matplotlib.pyplot as plt

def load_tsp(filename):
    points = {}
    with open(filename, "r") as f:
        start_reading = False
        for line in f:
            line = line.strip()
            if line == "NODE_COORD_SECTION":
                # print("start")
                start_reading = True
                continue
            if line == "EOF":
                # print("eof")
                break
            if start_reading:
                parts = line.split()
                if len(parts) < 3:
                    continue
                parts[0]=int(parts[0])
                points[parts[0]] = (int(parts[1]), int(parts[2]))
    return points


def load_results(filename):
    with open(filename, "r") as f:
        lines = f.readlines()
    indexes1 = lines[0].split()
    indexes2 = lines[1].split()
    return indexes1, indexes2

from math import sqrt
import numpy as np
def dist(p1,p2):
    return np.floor(sqrt((p1[0]-p2[0])**2+(p1[1]-p2[1])**2)+0.5)

def cycle_length(indexes, points):
    length = 0
    for i in range(len(indexes)):
        current_point = points[int(indexes[i])]
        next_point = points[int(indexes[(i + 1) % len(indexes)])]
        length += dist(current_point, next_point)
    return length

points = load_tsp("kroA200.tsp")
indexes1, indexes2 = load_results("trw.txt")

cycle1_length = cycle_length(indexes1, points)
cycle2_length = cycle_length(indexes2, points)
print(cycle1_length+cycle2_length)

fig, axs = plt.subplots(1, 2, figsize=(20, 10))

# print(points[0])
# print(points[183][0])

for idx in indexes1:
    axs[0].scatter(points[int(idx)][0],points[int(idx)][1], color="red")
for idx in indexes2:
    axs[0].scatter(points[int(idx)][0],points[int(idx)][1], color="blue")

for i in range(len(indexes1)):
    axs[0].plot([points[int(indexes1[i])][0], points[int(indexes1[(i+1)%len(indexes2)])][0]],
             [points[int(indexes1[i])][1], points[int(indexes1[(i+1)%len(indexes2)])][1]], color="red")

for i in range(len(indexes2)):
    axs[0].plot([points[int(indexes2[i])][0], points[int(indexes2[(i+1)%len(indexes2)])][0]],
             [points[int(indexes2[i])][1], points[int(indexes2[(i+1)%len(indexes2)])][1]], color="blue")
axs[0].set_title("kroA200.tsp")
axs[1].set_title("kroB200.tsp")

points = load_tsp("kroB200.tsp")
indexes1, indexes2 = load_results("trw2.txt")

cycle1_length = cycle_length(indexes1, points)
cycle2_length = cycle_length(indexes2, points)
print(cycle1_length+cycle2_length)

# print(points[0])
# print(points[183][0])

for idx in indexes1:
    axs[1].scatter(points[int(idx)][0],points[int(idx)][1], color="red")
for idx in indexes2:
    axs[1].scatter(points[int(idx)][0],points[int(idx)][1], color="blue")

for i in range(len(indexes1)):
    axs[1].plot([points[int(indexes1[i])][0], points[int(indexes1[(i+1)%len(indexes2)])][0]],
             [points[int(indexes1[i])][1], points[int(indexes1[(i+1)%len(indexes2)])][1]], color="red")

for i in range(len(indexes2)):
    axs[1].plot([points[int(indexes2[i])][0], points[int(indexes2[(i+1)%len(indexes2)])][0]],
             [points[int(indexes2[i])][1], points[int(indexes2[(i+1)%len(indexes2)])][1]], color="blue")

fig.suptitle("Candidates, k=10")
plt.show()