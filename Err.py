import pandas as pd
import matplotlib.pyplot as plt

gamm = 1.4

global file_path
global data
global time, Ratio_rho, Ratio_v, Ratio_P



def ReadingFile(name):
    global file_path
    global data
    global time, Ratio_rho, Ratio_v, Ratio_P

    file_path = "CSV_files/" + "Analis_"+ name + ".csv"  # путь к вашему файлу
    data = pd.read_csv(file_path)
    time = data['time']
    Ratio_rho = data['Ratio_rho']
    Ratio_v = data['Ratio_v']
    Ratio_P = data['Ratio_P']

colors = ['blue', 'red','purple', '-.k']
global n
n = 0

def FillAxes():
    global file_path
    global data
    global time, Ratio_rho, Ratio_v, Ratio_P
    global n
    global alpha

    axs[0].plot(time, Ratio_rho, colors[n], alpha=alpha)
    axs[0].set_xlabel('time')
    axs[0].set_ylabel('Ratio_rho')
    
    axs[1].plot(time,Ratio_v, colors[n], alpha=alpha)
    axs[1].set_xlabel('time')
    axs[1].set_ylabel('Ratio_v')

    axs[2].plot(time, Ratio_P, colors[n], alpha=alpha)
    axs[2].set_xlabel('time')
    axs[2].set_ylabel('Ratio_P')


plt.style.use('seaborn-v0_8')
fig, axs = plt.subplots(1, 3, figsize=(9, 7), sharex=True)
global alpha
alpha = 1

ReadingFile('Godunov')
FillAxes()

"""
ReadingFile('Riemann')
FillAxes()

ReadingFile('Kolgan')
FillAxes()

#ReadingFile('Rodionov')
#FillAxes()

ReadingFile('ENO')
FillAxes()

#ReadingFile('WENO')
#FillAxes()

#ReadingFile('acoust')
#FillAxes()

alpha = 0.6   

"""
# Единственная легенда вне графиков
#labels = ['Godunov', 'Kolgan', 'ENO', 'Analytical']
#fig.legend(labels, loc='center left', bbox_to_anchor=(0.9, 0.5), fontsize=12)

for ax in axs.flatten():
    ax.grid(True, alpha=0.5)

plt.tight_layout(rect=[0, 0, 0.85, 1])
plt.savefig(f'output/picture.png', bbox_inches='tight', dpi=300, transparent=False)
plt.show()

