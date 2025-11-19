# filter_bus_proper.py
import struct
import sys

target_ecs = {35213}

with open('output.unfiltered.bus', 'rb') as infile:
    # Read BUS header
    magic = infile.read(4)
    version = struct.unpack('<I', infile.read(4))[0]
    bclen = struct.unpack('<I', infile.read(4))[0]
    umilen = struct.unpack('<I', infile.read(4))[0]
    
    print(f"BUS version: {version}, BC length: {bclen}, UMI length: {umilen}")
    
    with open('filtered_GQX67_05945.bus', 'wb') as outfile:
        # Write header
        outfile.write(magic)
        outfile.write(struct.pack('<I', version))
        outfile.write(struct.pack('<I', bclen))
        outfile.write(struct.pack('<I', umilen))
        
        kept = 0
        total = 0
        
        # Read records (each is 32 bytes)
        while True:
            record = infile.read(32)
            if len(record) < 32:
                break
            
            total += 1
            
            # Extract EC (bytes 16-20 in the record)
            ec = struct.unpack('<I', record[16:20])[0]
            
            if ec in target_ecs:
                outfile.write(record)
                kept += 1
        
        print(f"Kept {kept} out of {total} records")
