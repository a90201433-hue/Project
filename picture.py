import pandas as pd
import matplotlib.pyplot as plt

gamm = 1.4

global file_path
global data
global x, rho, u, P, e



def ReadingFile(name):
    global file_path
    global data
    global x, rho, u, P, e

    file_path = "CSV_files/" + name + ".csv"  # путь к вашему файлу
    data = pd.read_csv(file_path)
    x = data['x']
    rho = data['rho']
    u = data['u']
    P = data['P']
    e = data['e']


colors = ['purple', '-.k']
global n
n = 0

def FillAxes():
    global file_path
    global data
    global x, rho, u, P, e
    global n
    global alpha

    axs[0, 0].plot(x, rho, colors[n], alpha=alpha)
    axs[0, 0].set_xlabel('x')
    axs[0, 0].set_ylabel(r'$\rho$')
    
    axs[0, 1].plot(x, u, colors[n], alpha=alpha)
    axs[0, 1].set_xlabel('x')
    axs[0, 1].set_ylabel('u')

    axs[1, 0].plot(x, P, colors[n], alpha=alpha)
    axs[1, 0].set_xlabel('x')
    axs[1, 0].set_ylabel('P')

    axs[1, 1].plot(x, e, colors[n], alpha=alpha)
    axs[1, 1].set_xlabel('x')
    axs[1, 1].set_ylabel(r'$\varepsilon$')
    n += 1


plt.style.use('seaborn-v0_8')
fig, axs = plt.subplots(2, 2, figsize=(9, 7), sharex=True)
global alpha
alpha = 1

#ReadingFile('Godunov')
#FillAxes()

#ReadingFile('Kolgan')
#FillAxes()

#ReadingFile('Rodionov')
#FillAxes()

#ReadingFile('ENO')
#FillAxes()

ReadingFile('WENO')
FillAxes()

#ReadingFile('acoust')
#FillAxes()

alpha = 0.6   
ReadingFile('Riemann')
FillAxes()

# Единственная легенда вне графиков
labels = ['ENO', 'Analytical']
fig.legend(labels, loc='center left', bbox_to_anchor=(0.9, 0.5), fontsize=12)

for ax in axs.flatten():
    ax.grid(True, alpha=0.5)

plt.tight_layout(rect=[0, 0, 0.85, 1])
plt.savefig(f'output/picture.png', bbox_inches='tight', dpi=300, transparent=False)
plt.show()

