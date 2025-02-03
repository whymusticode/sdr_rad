
# ssh root@192.168.2.1
# password: analog 

# from pluto.pluto_sdr import PlutoSdr
# sdr = PlutoSdr()

# watch -t SoapySDRUtil --find

# SoapySDRUtil --probe="driver=rtlsdr"
# rtl_test -d 0

# ip route show | grep 192.168.2


# [INFO] Auto setting Buffer Size: 65536
# [INFO] Set MTU Size: 65536

# /etc/modprobe.d/blacklist-rtl.conf
# blacklist dvb_usb_rtl28xxu
# blacklist rtl2832
# blacklist rtl2830


# # unload/reload driver 
# sudo rmmod rtl2832_sdr
# sudo rmmod dvb_usb_rtl28xxu
# sudo modprobe rtl2832_sdr




from gnuradio import gr, soapy
import time
import usb.core
import usb.util
from multiprocessing import Process, Queue
from gnuradio import gr, soapy, blocks
import time
import os
os.chdir('/home/bill/BACKUP/CODE/SDR/')
from scipy.fft import fft, ifft
import matplotlib.pyplot as plt
import importlib
import util
import numpy as np
importlib.reload(util)


def test_rtl(device_path):
    tb = gr.top_block()
    dev_args = f'driver=rtlsdr,device={device_path}'
    
    try:
        sdr = soapy.source(dev_args, "fc32", 1, '',
                          'bufflen=16384', [''], [''])
        sdr.set_sample_rate(0, 2.048e6)
        sdr.set_frequency(0, 99.5e6)
        sdr.set_gain(0, 'TUNER', 20.0)
        print(f"Successfully initialized device {device_path}")
        tb.start()
        time.sleep(1)
        tb.stop()
        return True
    except Exception as e:
        print(f"Error with device {device_path}: {e}")
        return False
    finally:
        tb = None
# Test each device individually using USB paths
print("Testing first device...")
test_rtl('0')
time.sleep(2)
print("\nTesting second device...")
test_rtl('1')






def sdr_process(queue, frequency, sample_rate, device_path):
    try:
        # Reset USB device first
        reset_usb_device()
        time.sleep(1)
        
        tb = gr.top_block()
        dev_args = f'driver=rtlsdr,device={device_path}'
        
        # Add initialization status message
        queue.put(f"Initializing device {device_path}")
        
        sdr = soapy.source(dev_args, "fc32", 1, '',
                          'bufflen=16384', [''], [''])
        
        sdr.set_sample_rate(0, sample_rate)
        sdr.set_frequency(0, frequency)
        sdr.set_gain(0, 'TUNER', 20.0)
        
        sink = blocks.vector_sink_c()
        tb.connect(sdr, sink)
        
        # Add streaming status message
        queue.put(f"Starting flowgraph for device {device_path}")
        
        tb.start()
        
        # Wait for actual data
        timeout = 0
        while timeout < 50:  # 5 second timeout
            time.sleep(0.1)
            data = sink.data()
            if len(data) > 0:
                queue.put(data[:1024])
                break
            timeout += 1
            
        if timeout >= 50:
            queue.put(f"Timeout waiting for data from device {device_path}")
            
    except Exception as e:
        queue.put(f"Error in device {device_path}: {str(e)}")
        raise e
    finally:
        if 'tb' in locals():
            tb.stop()
            tb.wait()

# Usage
q1, q2 = Queue(), Queue()
frequency = 99.5e6
sample_rate = 2.048e6

# Start first device
p1 = Process(target=sdr_process, args=(q1, frequency, sample_rate, '0'))
p1.start()

# Get initialization messages and data
while True:
    msg = q1.get()
    print(f"Device 0:", msg)
    if isinstance(msg, (list, tuple, np.ndarray)):
        samples1 = msg
        break

time.sleep(2)

# Start second device
p2 = Process(target=sdr_process, args=(q2, frequency, sample_rate, '1'))
p2.start()

# Get initialization messages and data with timeout
timeout_start = time.time()
while time.time() - timeout_start < 10:  # 10 second timeout
    try:
        msg = q2.get(timeout=1)
        print(f"Device 1:", msg)
        if isinstance(msg, (list, tuple, np.ndarray)):
            samples2 = msg
            break
    except queue.Empty:
        print("Timeout waiting for device 1")
        break





def reset_usb_device(vendor_id=0x0bda, product_id=0x2838):
    """Reset USB device before attempting to use it"""
    dev = usb.core.find(idVendor=vendor_id, idProduct=product_id)
    if dev is not None:
        try:
            dev.reset()
        except:
            pass
        usb.util.dispose_resources(dev)

def sdr_process(queue, frequency, sample_rate, device_path):
    try:
        # Reset USB device first
        reset_usb_device()
        time.sleep(1)  # Give USB time to re-enumerate
        
        tb = gr.top_block()
        dev_args = f'driver=rtlsdr,device={device_path}'
        
        sdr = soapy.source(dev_args, "fc32", 1, '',
                          'bufflen=16384', [''], [''])
        
        sdr.set_sample_rate(0, sample_rate)
        sdr.set_frequency(0, frequency)
        sdr.set_gain(0, 'TUNER', 20.0)
        
        sink = blocks.vector_sink_c()
        tb.connect(sdr, sink)
        
        tb.start()
        
        while True:
            time.sleep(0.1)
            data = sink.data()
            sink.reset()
            if len(data) > 0:
                queue.put(data[:1024])
                
    except Exception as e:
        queue.put(f"Error in device {device_path}: {str(e)}")
        raise e
    finally:
        if 'tb' in locals():
            tb.stop()
            tb.wait()

# Queue setup
q1, q2 = Queue(), Queue()

frequency = 99.5e6
sample_rate = 2.048e6

# Reset all devices first
reset_usb_device()
time.sleep(2)  # Give USB time to settle

p1 = Process(target=sdr_process, args=(q1, frequency, sample_rate, '0'))
p1.start()
samples1 = q1.get()

time.sleep(2)  # Longer delay between processes

p2 = Process(target=sdr_process, args=(q2, frequency, sample_rate, '1'))
p2.start()
samples2 = q2.get()




plt.ion()
fig, ax = plt.subplots()
line, = ax.plot([], [], lw=2)
ax.set_xlim(0, 1000)
ax.set_ylim(-20, 20)

ct = 0
while ct < 32:
    # if not q1.empty() and not q2.empty():
        samples1 = q1.get()
        samples2 = q2.get()
        ct += 1
        samples1 = np.array(samples1).astype(np.complex64)
        samples2 = np.array(samples2).astype(np.complex64)
        dbDat = util.db(fft(samples1))

        line.set_data(range(len(dbDat)), dbDat)
        plt.pause(0.1)



p1.terminate()
p2.terminate()
p1.join()
p2.join()












def sdr_process(queue, frequency, sample_rate, device_path):
    tb = gr.top_block()
    
    # Specify device by index
    # dev_args = f'driver=rtlsdr,rtl={device_index}'
    dev_args = f'driver=rtlsdr,device={device_path}'
    
    sdr = soapy.source(dev_args, "fc32", 1, '',
                       'bufflen=16384', [''], [''])
    
    # Configure SDR
    sdr.set_sample_rate(0, sample_rate)
    sdr.set_frequency(0, frequency)
    sdr.set_gain(0, 'TUNER', 20.0)
    
    # Create sink to capture data
    sink = blocks.vector_sink_c()
    tb.connect(sdr, sink)
    
    # Start flowgraph
    tb.start()
    
    while True:
        time.sleep(0.1)
        data = sink.data()
        sink.reset()
        if len(data) > 0:
            queue.put(data[:1024])



from gnuradio import gr, soapy
tb = gr.top_block()
# sdr = soapy.source('driver=rtlsdr,rtl=0', "fc32", 1, '', 'bufflen=16384', [''], [''])
# tb.start()






from multiprocessing import Process, Queue
from gnuradio import gr, soapy, blocks
import numpy as np
import time




# Queue setup
q1, q2 = Queue(), Queue()

frequency = 99.5e6
sample_rate = 2.048e6

# Create processes using different device indices
p1 = Process(target=sdr_process, args=(q1, frequency, sample_rate, '0'))  # First device
p2 = Process(target=sdr_process, args=(q2, frequency, sample_rate, '1'))  # Second device

p1.start()
# time.sleep(1)  # Add delay between starting processes
p2.start()






from multiprocessing import Process, Queue
from gnuradio import gr, soapy, blocks
import time

def sdr_process(queue, frequency, sample_rate, device_path):
    try:
        tb = gr.top_block()
        dev_args = f'driver=rtlsdr,device={device_path}'
        
        sdr = soapy.source(dev_args, "fc32", 1, '',
                          'bufflen=16384', [''], [''])
        
        sdr.set_sample_rate(0, sample_rate)
        sdr.set_frequency(0, frequency)
        sdr.set_gain(0, 'TUNER', 20.0)
        
        sink = blocks.vector_sink_c()
        tb.connect(sdr, sink)
        
        # Add a small delay to stagger device initialization
        time.sleep(float(device_path) * 0.5)
        
        tb.start()
        
        while True:
            time.sleep(0.1)
            data = sink.data()
            sink.reset()
            if len(data) > 0:
                queue.put(data[:1024])
                
    except Exception as e:
        queue.put(f"Error in device {device_path}: {str(e)}")
        raise e
    finally:
        if 'tb' in locals():
            tb.stop()
            tb.wait()

# Queue setup
q1, q2 = Queue(), Queue()

frequency = 99.5e6
sample_rate = 2.048e6

# Create and start processes with staggered timing
p1 = Process(target=sdr_process, args=(q1, frequency, sample_rate, '0'))
p2 = Process(target=sdr_process, args=(q2, frequency, sample_rate, '1'))

p1.start()
time.sleep(1)  # Add significant delay between starting processes
p2.start()





def sdr_process(queue, frequency, sample_rate, device_serial=''):
    tb = gr.top_block()
    
    # Specify the device by serial number
    dev_args = f'driver=rtlsdr,serial={device_serial}' if device_serial else 'driver=rtlsdr'
    
    sdr = soapy.source(dev_args, "fc32", 1, '',
                       'bufflen=16384', [''], [''])
    
    # Configure SDR
    sdr.set_sample_rate(0, sample_rate)
    sdr.set_frequency(0, frequency)
    sdr.set_gain(0, 'TUNER', 20.0)
    
    # Create sink to capture data
    sink = blocks.vector_sink_c()
    tb.connect(sdr, sink)
    
    # Start flowgraph
    tb.start()
    
    while True:
        time.sleep(0.1)
        data = sink.data()
        sink.reset()
        if len(data) > 0:
            queue.put(data[:1024])

# First, let's get the serial numbers of our RTL-SDR devices
def get_rtl_serials():
    import subprocess
    try:
        result = subprocess.run(['rtl_test'], capture_output=True, text=True)
        # Parse the output to get device serials
        # You might need to adjust this parsing based on your rtl_test output
        return ['00000001', '00000002']  # Replace with actual serials
    except Exception as e:
        print(f"Error getting RTL-SDR serials: {e}")
        return []

# Get device serials
rtl_serials = get_rtl_serials()
if len(rtl_serials) < 2:
    raise RuntimeError("Not enough RTL-SDR devices found")

# Queue setup
q1, q2 = Queue(), Queue()

frequency = 99.5e6
sample_rate = 2.048e6

# Create and start processes with specific device serials
p1 = Process(target=sdr_process, args=(q1, frequency, sample_rate, rtl_serials[0]))
p2 = Process(target=sdr_process, args=(q2, frequency, sample_rate, rtl_serials[1]))

try:
    p1.start()
    time.sleep(1)  # Add delay between starting processes
    p2.start()

    plt.ion()
    fig, ax = plt.subplots()
    line, = ax.plot([], [], lw=2)
    ax.set_xlim(0, 1000)
    ax.set_ylim(-20, 20)

    ct = 0
    while ct < 32:
        if not q1.empty() and not q2.empty():
            samples1 = q1.get()
            samples2 = q2.get()
            ct += 1
            samples1 = np.array(samples1).astype(np.complex64)
            samples2 = np.array(samples2).astype(np.complex64)
            dbDat = util.db(fft(samples1))
            line.set_data(range(len(dbDat)), dbDat)
            plt.pause(0.1)

finally:
    # Cleanup
    p1.terminate()
    p2.terminate()
    p1.join()
    p2.join()


    

import os
os.chdir('/home/bill/BACKUP/CODE/SDR/')

from multiprocessing import Process, Queue
from gnuradio import gr, soapy, blocks
import numpy as np
import time
from scipy.fft import fft, ifft
import matplotlib.pyplot as plt
import importlib
import util

def sdr_process(queue, frequency, sample_rate, uri=''):
    tb = gr.top_block()
    
    # Modified for RTL-SDR
    sdr = soapy.source('driver=rtlsdr', "fc32", 1, uri, 
                       'bufflen=16384', [''], [''])
    
    # Configure SDR - adjusted for RTL-SDR parameters
    sdr.set_sample_rate(0, sample_rate)
    sdr.set_frequency(0, frequency)
    sdr.set_gain(0, 'TUNER', 20.0)  # RTL-SDR uses 'TUNER' gain
    
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

# Queue setup
q1, q2 = Queue(), Queue()

# Modified frequency for RTL-SDR range (87.5 MHz - 108 MHz for FM radio)
frequency = 99.5e6  # Example: 99.5 MHz FM radio
sample_rate = 2.048e6  # RTL-SDR supported rate

# Create and start processes
p1 = Process(target=sdr_process, args=(q1, frequency, sample_rate, ''))
p2 = Process(target=sdr_process, args=(q2, frequency, sample_rate, ''))

p1.start()
p2.start()

# Plotting setup
plt.ion()  # Turn on interactive mode
fig, ax = plt.subplots()
line, = ax.plot([], [], lw=2)
ax.set_xlim(0, 1000)
ax.set_ylim(-20, 20)

ct = 0
while ct < 32:
    if not q1.empty() and not q2.empty():
        samples1 = q1.get()
        samples2 = q2.get()
        ct += 1
        samples1 = np.array(samples1).astype(np.complex64)
        samples2 = np.array(samples2).astype(np.complex64)
        dbDat = util.db(fft(samples1))
        line.set_data(range(len(dbDat)), dbDat)
        plt.pause(0.1)

# Clean up
p1.terminate()
p2.terminate()
p1.join()
p2.join()




import os
os.chdir('/home/bill/BACKUP/CODE/SDR/')

from multiprocessing import Process, Queue
from gnuradio import gr, soapy, blocks
import numpy as np
import time
from scipy.fft import fft, ifft
import matplotlib.pyplot as plt


import importlib

import util
importlib.reload(util)
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
        # util.fftConv(samples1,samples2,flag='full')
        dbDat = util.db(fft(samples1))
        line.set_data(range(len(dbDat)), dbDat)  # Update the line data
        plt.pause(0.1)       # Pause to update the plot
        

# Nstep = 256

# import plotly.graph_objects as go
# from plotly.subplots import make_subplots

# # Create figure once
# fig = go.Figure()
# fig.add_trace(go.Scatter(y=db(util.randn_iq(100))))
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



ambig = util.fftConv(data[:,0,0],data[:,0,1],'full')

thing = fft(data[:,0,0])
util.quickPlot(thing)


# from gnuradio import soapy
# stream_args = ''
# tune_args = ['']
# settings = ['']
# sdr1 = soapy.source('driver=plutosdr', "fc32", 1, '', stream_args, tune_args, settings)






