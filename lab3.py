import numpy as np
from PIL import Image
from matplotlib import pyplot as plt
import komm


def close_event():
    plt.close()
    plt.close()

def qam_bases(n):
    n_val = np.sqrt(n)/2
    if n_val != int(n_val):
        raise AssertionError("Not 2^2m")
    n_val = int(n_val)
    vals = np.linspace(1,2*n_val-1,n_val)
    print(vals)
    total_e = 0
    for q in vals:
        for i in vals:
            print(i,q)
            total_e += q**2+i**2
    e_sym = total_e/len(vals)**2
    print(e_sym)
    vals = vals/np.sqrt(e_sym)
    total_e = 0
    for q in vals:
        for i in vals:
            print(i,q)
            total_e += q**2+i**2
    e_sym = total_e/len(vals)**2
    print(e_sym)
    return vals[0]

def calc_pad(length,sym_size):
    return (sym_size-(length%sym_size))%sym_size

m = 64
n = qam_bases(m)
print(n)
modulations = [komm.PSKModulation(2),komm.PSKModulation(4,phase_offset=np.pi/4)] + [(komm.QAModulation(16*2**(2*n),qam_bases(16*2**(2*n)))) for n in range(3)]
print(modulations)
tx_im = Image.open("imsend-10.pgm")
Npixels = tx_im.size[1]*tx_im.size[0]
plt.figure()
plt.imshow(np.array(tx_im),cmap="gray",vmin=0,vmax=255)
plt.show()
tx_bin_raw = np.unpackbits(np.array(tx_im))
n_bits = len(tx_bin_raw)


ber_ray = []
ber_db = range(-10,20)
for mod in modulations:
    pad = calc_pad(n_bits,mod.bits_per_symbol)
    tx_bin = np.append(tx_bin_raw,np.zeros(pad))
    for i in ber_db:
        print(mod.energy_per_symbol)
        awgn = komm.AWGNChannel(snr=10**(i/10))
        tx_data = mod.modulate(tx_bin)
        rx_data = awgn(tx_data)
        rx_bin = mod.demodulate(rx_data)
        if pad == 0:
            rx_bin_raw = rx_bin
        else:
            rx_bin_raw = rx_bin[:-pad]
        print(tx_bin_raw,rx_bin_raw)
        cmp = tx_bin_raw != rx_bin_raw
        errors = sum(cmp)
        bit_error_rate = errors/Npixels/8
        ber_ray.append(bit_error_rate)
        rx_im = np.packbits(rx_bin_raw).reshape(tx_im.size[1],tx_im.size[0])
        fig = plt.figure()
        timer = fig.canvas.new_timer(interval = 700) #creating a timer object and setting an interval of 3000 milliseconds
        timer.add_callback(close_event)
        plt.imshow(np.array(rx_im),cmap="gray",vmin=0,vmax=255)
        timer.start()
        plt.figure()
        plt.axes().set_aspect("equal")
        plt.scatter(rx_data[:100000].real,rx_data[:100000].imag,s=1,marker=".")
        plt.show()
    timer.stop()
    plt.plot(ber_db,ber_ray)
    plt.yscale("log")
    plt.show()
    timer.start()



