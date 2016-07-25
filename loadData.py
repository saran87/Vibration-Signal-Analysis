import json
import math
import numpy as np


def splitAndConvert(value, splitter=' '):
    return map(float, value.split(splitter))


def getRMSData(x, y, z):
    return [math.sqrt((x[i] * x[i]) + (y[i] * y[i]) + (z[i] * z[i])) for i in xrange(0, 2048)]


def loadEvents():
    events = []
    with open('vibration.json') as data_file:
        for line in data_file:
            event = json.loads(line)
            event['x'] = splitAndConvert(event['value']['x']);
            event['y'] = splitAndConvert(event['value']['y']);
            event['z'] = splitAndConvert(event['value']['z']);
            event['rms'] = getRMSData(event['x'], event['y'], event['z'])
            del event['value'];
            events.append(event)
    return events


def processPSD(events, axis):
    psds = []
    for event in events:
        signal = np.array(event[axis], dtype=float)
        fourier = np.fft.fft(signal * np.hanning(2048))
        psd = (fourier.real * fourier.real) / (1600 * 2048)
        psds.append(psd)
    return psds


def consolidatePSD(psds):
    dArray = np.array(psds)
    print 'Data type                :', dArray.dtype
    print 'Total number of elements :', dArray.size
    print 'Number of dimensions     :', dArray.ndim
    print 'Shape (dimensionality)   :', dArray.shape
    print 'Memory used (in bytes)   :', dArray.nbytes
    vibrationTable = []
    for i in xrange(0, 2048):
        vibrationTable.append(np.sum(dArray[:, i]) / len(psds))

    return vibrationTable


def processData():
    events = loadEvents()
    table = {}
    table['x'] = consolidatePSD(processPSD(events, 'x'))
    table['y'] = consolidatePSD(processPSD(events, 'y'))
    table['z'] = consolidatePSD(processPSD(events, 'z'))
    return table
