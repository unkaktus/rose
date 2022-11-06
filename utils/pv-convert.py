#!/usr/bin/env python
import paraview.simple as pv
import argparse
import os

parser = argparse.ArgumentParser('rose')
parser.add_argument("-i", type=str, help="Input filename")
parser.add_argument("-o", type=str, help="Output filename")

args = parser.parse_args()

reader = pv.OpenDataFile(args.i)
writer = pv.CreateWriter(args.o, reader)
writer.WriteAllTimeSteps = 1
writer.UpdatePipeline()