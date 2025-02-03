
# ssh root@192.168.2.1
# password: analog 

# from pluto.pluto_sdr import PlutoSdr
# sdr = PlutoSdr()

# watch -t SoapySDRUtil --find
# ip route show | grep 192.168.2


# [INFO] Auto setting Buffer Size: 65536
# [INFO] Set MTU Size: 65536

import os
os.chdir('/home/bill/BACKUP/CODE/SDR/')

from multiprocessing import Process, Queue
from gnuradio import gr, soapy, blocks
import numpy as np
import time
from scipy.fft import fft, ifft
import matplotlib.pyplot as plt


import importlib

import eca2
importlib.reload(eca2)
# from IPython import get_ipython
# get_ipython().run_line_magic('load_ext', 'autoreload')
# get_ipython().run_line_magic('autoreload', '2')



def sdr_process(queue, frequency, sample_rate, uri=''):
    tb = gr.top_block()
    sdr = soapy.source('driver=plutosdr', "fc32", 1, uri, '', [''], [''])
    
    # Configure SDR
    sdr.set_sample_rate(0, sample_rate)
    sdr.set_frequency(0, frequency)
    sdr.set_gain(0, 20.0)
    
    # Create sink to capture data
    sink = blocks.vector_sink_c()
    tb.connect(sdr, sink)
    
    # Start flowgraph
    tb.start()
    
    while True:
        time.sleep(0.1)  # Collect for 100ms
        data = sink.data()
        sink.reset()
        if len(data) > 0:
            queue.put(data[:1024])  # Send 1024 samples


q1, q2 = Queue(), Queue()

frequency = 2.437e9   # 99.5e6
sample_rate = 2.048e6
# Create and start processes
p1 = Process(target=sdr_process, args=(q1,frequency, sample_rate, ''))
p2 = Process(target=sdr_process, args=(q2,frequency, sample_rate, ''))

p1.start()
p2.start()




plt.ion()  # Turn on interactive mode
fig, ax = plt.subplots()
line, = ax.plot([], [], lw=2)
ax.set_xlim(0, 1000)
ax.set_ylim(-20, 20)

ct = 0
while ct<32:
    if not q1.empty() and not q2.empty():
        samples1 = q1.get()
        samples2 = q2.get()
        ct += 1
        samples1 = np.array(samples1).astype(np.complex64)
        samples2 = np.array(samples2).astype(np.complex64)
        # eca2.fftConv(samples1,samples2,flag='full')
        dbDat = eca2.db(fft(samples1))
        line.set_data(range(len(dbDat)), dbDat)  # Update the line data
        plt.pause(0.1)       # Pause to update the plot
        

# Nstep = 256

# import plotly.graph_objects as go
# from plotly.subplots import make_subplots

# # Create figure once
# fig = go.Figure()
# fig.add_trace(go.Scatter(y=db(eca2.randn_iq(100))))
# fig.update_layout(template='plotly_dark', hovermode=False)
# fig.show()

# # Initialize figure
# fig = go.Figure()
# fig.add_trace(go.Scatter(y=np.zeros(100), mode="lines", name="Channel 1"))

# # Display the initial figure (if in a Jupyter Notebook or Dash app)
# fig.show()


# # Set axis limits


# Simulate dynamic data
for _ in range(32):
    x = np.arange(100)
    y = 10 * np.log10(np.abs(np.random.random(100)) + 1e-9)
    line.set_data(x, y)  # Update the line data
    plt.pause(0.1)       # Pause to update the plot




import dash
from dash import dcc, html
from dash.dependencies import Input, Output
import numpy as np
import plotly.graph_objects as go

app = dash.Dash(__name__)
app.layout = html.Div([
    dcc.Graph(id="live-plot"),
    dcc.Interval(id="interval", interval=100, n_intervals=0)
])

def db(samples):
    return 10 * np.log10(np.abs(samples) + 1e-9)

@app.callback(
    Output("live-plot", "figure"),
    [Input("interval", "n_intervals")]
)
def update_graph(n):
    x = np.arange(100)
    y = db(np.random.random(100))
    fig = go.Figure(go.Scatter(x=x, y=y, mode="lines"))
    return fig

app.run_server(debug=True)




# while True:
#     if not q1.empty() and not q2.empty():
#         samples1 = q1.get()
#         samples2 = q2.get()
#         print()

def quickPlot(thing):
    dbDat = np.log10(np.real(thing*np.conj(thing)))*10
    fig = go.Figure(data=go.Scatter(y=dbDat))
    fig.update_layout(template='plotly_dark',hovermode=False)
    fig.show()

ct = 0
Nstep = 128
data = np.zeros((1024,2,Nstep)).astype(np.complex64)
while ct < Nstep:
    if not q1.empty() and not q2.empty():
        samples1 = q1.get()
        samples2 = q2.get()
        data[:,0,ct] = samples1
        data[:,1,ct] = samples2
        ct += 1
        print(ct)
        quickPlot(fft(samples1))
        # print(f"Got samples - SDR1: {len(samples1)}, SDR2: {len(samples2)}")


import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np

# Create figure and axis
fig, ax = plt.subplots()
line, = ax.plot([], [])

# Initialize plot limits
ax.set_xlim(0, 100)
ax.set_ylim(-1, 1)

def update(frame):
    # Get your new data here
    x = np.linspace(0, 100, 100)
    y = np.sin(x/10 + frame/10)
    line.set_data(x, y)
    return line,

ani = FuncAnimation(fig, update, frames=None, interval=50)
plt.show()



ambig = eca2.fftConv(data[:,0,0],data[:,0,1],'full')

thing = fft(data[:,0,0])
eca2.quickPlot(thing)


# from gnuradio import soapy
# stream_args = ''
# tune_args = ['']
# settings = ['']
# sdr1 = soapy.source('driver=plutosdr', "fc32", 1, '', stream_args, tune_args, settings)






