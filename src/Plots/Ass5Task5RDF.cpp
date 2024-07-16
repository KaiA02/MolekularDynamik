//code für diffusion output
//import matplotlib.pyplot as plt
//
//# Read the file (replace 'filename.txt' with your actual file name)
//        with open('rdfsupercooling.txt', 'r') as file:
//lines = file.readlines()
//
//# Extract intervals and local densities
//timesteps = []
//intervals = []
//local_densities = []
//for line in lines:
//parts = line.split()
//for i in range(1, len(parts)-1):
//intervals.append(float(parts[i]))
//local_densities.append(float(parts[i+1]))
//i+=1
//# Example data (replace with your actual data)
//tenth_lines = lines[::10]
//
//
//
//#Print the extracted data (optional)
//for line in tenth_lines:
//parts = line.split()
//timestep = int(parts[0])
//intervals = [float(val) for val in parts[1::2]]  # Extract every second value (intervals)
//local_densities = [float(val) for val in parts[2::2]]  # Extract every second value (local densities)
//plt.plot(intervals, local_densities, label=f'Timestep {timestep}')
//plt.xlabel('Intervals')
//plt.ylabel('Local Densities')
//plt.legend()
//plt.title('Local Densities vs. Intervals for Each Timestep')
//plt.grid(True)
//
//plt.show()

/**
 * @brief RDF Cooling
 * \image html RDFCooling.png width=800cm height=600cm
 */
void rdfCooling();

/**
 * @brief RDF SuperCooling
 * \image html RDFSuperCooling.png width=800cm height=600cm
 */
void rdfSuperCooling();