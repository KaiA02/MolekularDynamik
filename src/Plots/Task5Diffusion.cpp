//code für diffusion output
//import matplotlib.pyplot as plt
//
//# Funktion zum Einlesen der Daten aus einer .txt-Datei
//        def read_data(file_path):
//iterations = []
//values = []
//with open(file_path, 'r') as file:
//for line in file:
//parts = line.strip().split()
//try:
//iterations.append(int(parts[0]))
//values.append(float(parts[1]))
//except (IndexError, ValueError) as e:
//print(f"Fehler beim Verarbeiten der Zeile: {line.strip()} - {e}")
//return iterations, values
//
//# Daten aus den beiden Dateien einlesen
//        iterations1, values1 = read_data('diffusioncooling.txt')
//iterations2, values2 = read_data('diffusionsupercooling.txt')
//
//# Plot erstellen
//plt.figure(figsize=(10, 6))
//plt.plot(iterations1, values1, label='Cooling', color='blue', marker='o')
//plt.plot(iterations2, values2, label='SuperCooling', color='red', marker='x')
//
//# Plot-Details hinzufügen
//plt.xlabel('TimeStep')
//plt.title('Diffusion Var(t)')
//plt.legend()
//plt.grid(True)
//
//# Plot anzeigen
//plt.show()

/**
 * @brief Diffusion Var(t)
 * \image html Diffusion.png width=800cm height=600cm
 */
 void diffusion()