import numpy as np
import math
from decimal import Decimal
from tqdm import tqdm
def read_stl(output_vertices=False, filename="", output_facet_normal=True, int_vertex = np.asarray([0,0,0])):
    
    if isinstance(filename, np.ndarray):
        lines = filename
    else:
        with open(filename, "r") as f:

            lines = f.readlines()
    
    if output_vertices:
        vertices = []
        print("Calculating Vertices ... \n")
        
        for line in tqdm(lines):
            if isinstance(line, str):
                split_line = line.split()
                if split_line[0] == "vertex":
                    #vertex = float(split_line[1])
                    vertex = [float(split_line[1]),float(split_line[2]),float(split_line[3])]
                    vertices.append(vertex)
            else:
                vertices.append(line)
        
        return np.asarray(vertices)
    if output_facet_normal:
        vertices = []
        print("Calculating Vertices ... \n")

        for line in tqdm(lines):
            
            if isinstance(line, str):
                split_line = line.split()
                if split_line[0] == "facet":
                    facet_normal = [float(split_line[2]),float(split_line[3]),float(split_line[4])]
                elif split_line[0] == "vertex":
                    vertex = np.asarray([float(split_line[1]),float(split_line[2]),float(split_line[3])])
                    #vertex = [float(split_line[1]),float(split_line[2]),float(split_line[3])]
                    vertices.append((vertex, facet_normal))
            #else:
             #   vertices.append((line,facet_normal))

        for vertex, normal in vertices:

            if float("{:,.3f}".format(vertex[0])) == int_vertex[0] and float("{:,.3f}".format(vertex[1])) == int_vertex[1] and float("{:,.3f}".format(vertex[2])) == int_vertex[2]:
                return vertex, normal
        
        return np.asarray(vertices)
        
    return -1


if __name__ == "__main__":

    vertex, normal = read_stl(filename="493_STL_mRCA.stl", int_vertex=np.asarray([-58.434,-552.182,236.319]))
    
    print(f'X: {180*np.arccos(np.dot(np.asarray([1,0,0]),normal)) / 3.1415169}')
    print(f'Y: {180*np.arccos(np.dot(np.asarray([0,1,0]),normal)) / 3.1415169}')
    print(f'Z: {180*np.arccos(np.dot(np.asarray([0,0,1]),normal)) / 3.1415169}')
    print(vertex,normal)


