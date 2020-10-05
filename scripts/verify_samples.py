#!/usr/bin/env/python3

import csv
import os

with open('sampledata.csv', 'r', newline='') as infile:
  dialect = csv.Sniffer().sniff(infile.read())
  infile.seek(0)
  reader = csv.DictReader(infile, dialect=dialect)
  names = set()
  filenames = set()
  ids = set()
  ctr = 0
  first = True
  sadict = {'Row': [], 'Sample': [], 'Id': [], 'Name': []}
  print(sadict)

  for row in reader:
    print(row)
    if not 'Name' in row:
      raise RuntimeError('sampledata.csv needs a column for "Name"')
    if not 'Sample' in row:
      raise RuntimeError('sampledata.csv needs a column for "Sample"')
    if not 'BasespaceId' in row:
      raise RuntimeError('sampledata.csv needs a column for "BasespaceId"')
    if not row['Sample'] and not row['BasespaceId']:
      raise RuntimeError('Please provide either "Sample" or "BasespaceId"')
    sadict['Row'].append(ctr)
    ctr += 1
    if row['Name']:
      names.add(row['Name'])
      sadict['Name'].append(row['Name'])
    if row['Sample']:
     filenames.add(row['Sample'])
     sadict['Sample'].append(row['Sample'])
    if row['BasespaceId']:
      ids.add(row['BasespaceId'])
      sadict['Id'].append(row['BasespaceId'])

if ids and filenames:
  if not (len(filenames) == len(names) == len(ids) == ctr):
    raise RuntimeError('Columns "Name", "Sample" and "BasespaceId" must be unique')
  for i in range(ctr):
    print('BS', sadict['Sample'][i], sadict['Id'][i], sep='\t')

elif ids:
  if not (len(names) == len(ids) == ctr):
    raise RuntimeError('Columns "Name" and "BasespaceId" must be unique')
  for i in range(ctr):
    print('BS', sadict['Name'][i].strip() + '.fastq.gz', sadict['Id'][i], sep='\t')

elif filenames:
  if not (len(filenames) == len(names) == ctr):
    raise RuntimeError('Columns "Name" and "Sample" must be unique')
  for i in range(ctr):
    if not sadict['Sample'][i] in os.listdir('./rawdata'):
      raise RuntimeError(f"Sample {sadict['Sample'][i]} not found")
  for i in range(ctr):
    print('LOCAL', sadict['Sample'][i].strip(), sep='\t')
