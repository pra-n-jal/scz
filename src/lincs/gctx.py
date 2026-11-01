import h5py
import numpy as np
import pandas as pd
import sys
import os
from pathlib import Path
import argparse

def safe_read_gctx_top_lines(file_path, n_lines=10, max_memory_gb=2):
    """
    Safely read the top n lines of a GCTX file without loading everything into memory
    
    Args:
        file_path: Path to the GCTX file
        n_lines: Number of lines/rows to display (default 10)
        max_memory_gb: Maximum memory to use in GB (default 2)
    """
    
    # Check if file exists
    if not os.path.exists(file_path):
        print(f"Error: File not found at {file_path}")
        return
    
    # Check file size
    file_size_gb = os.path.getsize(file_path) / (1024**3)
    print(f"File size: {file_size_gb:.2f} GB")
    
    try:
        # Open the GCTX file safely
        print("Opening GCTX file...")
        with h5py.File(file_path, 'r') as gctx_file:
            # Print file structure first
            print("\n=== GCTX File Structure ===")
            def print_structure(name, obj):
                if isinstance(obj, h5py.Group):
                    print(f"Group: {name}")
                elif isinstance(obj, h5py.Dataset):
                    print(f"Dataset: {name}, Shape: {obj.shape}, Type: {obj.dtype}")
            
            gctx_file.visititems(print_structure)
            
            # Look for the main data matrix (usually in /0/DATA/0/matrix)
            matrix_path = "/0/DATA/0/matrix"
            if matrix_path in gctx_file:
                print(f"\n=== Found data matrix at {matrix_path} ===")
                matrix = gctx_file[matrix_path]
                print(f"Matrix shape: {matrix.shape}")
                print(f"Matrix dtype: {matrix.dtype}")
                
                # Calculate safe chunk size based on memory limit
                element_size = matrix.dtype.itemsize
                max_elements = int((max_memory_gb * 1024**3) / element_size)
                n_cols = matrix.shape[1] if len(matrix.shape) > 1 else 1
                safe_rows = min(n_lines, max_elements // n_cols, matrix.shape[0])
                
                print(f"\nReading top {safe_rows} rows...")
                
                # Read only the top n rows
                if len(matrix.shape) == 2:
                    top_data = matrix[:safe_rows, :]
                else:
                    top_data = matrix[:safe_rows]
                
                print(f"\n=== Top {safe_rows} rows of data matrix ===")
                print(f"Shape of extracted data: {top_data.shape}")
                
                # Display the data
                if len(top_data.shape) == 2:
                    # Show first few columns if too many
                    max_cols_display = 10
                    if top_data.shape[1] > max_cols_display:
                        print(f"Showing first {max_cols_display} columns out of {top_data.shape[1]}:")
                        display_data = top_data[:, :max_cols_display]
                    else:
                        display_data = top_data
                    
                    # Convert to pandas for better display
                    df = pd.DataFrame(display_data)
                    print(df)
                else:
                    print(top_data)
            
            # Look for column metadata (col_meta)
            col_meta_paths = ["/0/META/COL", "/0/META/col", "/META/COL", "/META/col"]
            col_meta_found = False
            
            for col_path in col_meta_paths:
                if col_path in gctx_file:
                    print(f"\n=== Found column metadata at {col_path} ===")
                    col_meta = gctx_file[col_path]
                    
                    if hasattr(col_meta, 'keys'):
                        print("Column metadata fields:")
                        for key in col_meta.keys():
                            field = col_meta[key]
                            if hasattr(field, 'shape'):
                                print(f"  {key}: shape {field.shape}, dtype {field.dtype}")
                                # Show first few values
                                if field.shape[0] > 0:
                                    sample_size = min(5, field.shape[0])
                                    sample_data = field[:sample_size]
                                    print(f"    Sample values: {sample_data}")
                    col_meta_found = True
                    break
            
            if not col_meta_found:
                print("\n=== Column metadata not found in expected locations ===")
                print("Available paths:")
                def list_all_paths(name, obj):
                    print(f"  {name}")
                gctx_file.visititems(list_all_paths)
            
            # Look for row metadata (row_meta)
            row_meta_paths = ["/0/META/ROW", "/0/META/row", "/META/ROW", "/META/row"]
            row_meta_found = False
            
            for row_path in row_meta_paths:
                if row_path in gctx_file:
                    print(f"\n=== Found row metadata at {row_path} ===")
                    row_meta = gctx_file[row_path]
                    
                    if hasattr(row_meta, 'keys'):
                        print("Row metadata fields:")
                        for key in row_meta.keys():
                            field = row_meta[key]
                            if hasattr(field, 'shape'):
                                print(f"  {key}: shape {field.shape}, dtype {field.dtype}")
                                # Show first few values
                                if field.shape[0] > 0:
                                    sample_size = min(5, field.shape[0])
                                    sample_data = field[:sample_size]
                                    print(f"    Sample values: {sample_data}")
                    row_meta_found = True
                    break
            
            if not row_meta_found:
                print("\n=== Row metadata not found in expected locations ===")
    
    except MemoryError:
        print("Error: Not enough memory to read the file. Try reducing n_lines or max_memory_gb.")
    except Exception as e:
        print(f"Error reading GCTX file: {str(e)}")
        print(f"Error type: {type(e).__name__}")

def extract_lincs_col_metadata(file_path, output_path="lincs_col_meta.tsv", max_memory_gb=2):
    """
    Extract LINCS column metadata from GCTX file and save as TSV
    
    Args:
        file_path: Path to the GCTX file
        output_path: Output TSV file path
        max_memory_gb: Maximum memory to use in GB
    """
    
    if not os.path.exists(file_path):
        print(f"Error: File not found at {file_path}")
        return False
    
    file_size_gb = os.path.getsize(file_path) / (1024**3)
    print(f"[INFO] File size: {file_size_gb:.2f} GB")
    print("[INFO] Reading GCTX column metadata...")
    
    try:
        with h5py.File(file_path, 'r') as gctx_file:
            print("[DEBUG] Exploring file structure for metadata...")
            
            # First, let's see the full structure
            all_paths = []
            def collect_paths(name, obj):
                all_paths.append(name)
                if isinstance(obj, h5py.Dataset):
                    print(f"[DEBUG] Dataset: {name}, Shape: {obj.shape}, Type: {obj.dtype}")
                elif isinstance(obj, h5py.Group):
                    print(f"[DEBUG] Group: {name}")
            
            gctx_file.visititems(collect_paths)
            
            # Look for column metadata in multiple possible locations
            col_meta_paths = [
                "/0/META/COL", "/0/META/col", "/META/COL", "/META/col",
                "/0/META/column", "/META/column", "/column", "/COL"
            ]
            
            col_meta_data = {}
            col_meta_found = False
            successful_path = None
            
            for col_path in col_meta_paths:
                if col_path in gctx_file:
                    print(f"[INFO] Found column metadata at {col_path}")
                    successful_path = col_path
                    col_meta = gctx_file[col_path]
                    
                    if hasattr(col_meta, 'keys'):
                        print(f"[INFO] Available metadata fields: {list(col_meta.keys())}")
                        
                        # Extract all metadata fields
                        for key in col_meta.keys():
                            field = col_meta[key]
                            if hasattr(field, 'shape') and len(field.shape) > 0:
                                try:
                                    print(f"[INFO] Processing field '{key}' with shape {field.shape}...")
                                    
                                    # Read the field data
                                    field_data = field[:]
                                    
                                    # Handle different data types
                                    if field_data.dtype.kind in ['S', 'U']:  # String data
                                        if field_data.dtype.kind == 'S':
                                            # Byte strings - decode to UTF-8
                                            field_data = [x.decode('utf-8') if isinstance(x, bytes) else str(x) for x in field_data]
                                        else:
                                            # Unicode strings
                                            field_data = [str(x) for x in field_data]
                                    elif field_data.dtype.kind in ['i', 'f']:  # Numeric data
                                        # Keep as numeric but convert to list for consistency
                                        field_data = field_data.tolist()
                                    else:
                                        # Other types - convert to string
                                        field_data = [str(x) for x in field_data]
                                    
                                    col_meta_data[key] = field_data
                                    print(f"[INFO] Extracted field '{key}': {len(field_data)} entries, sample: {field_data[:3] if len(field_data) > 0 else 'empty'}")
                                    
                                except Exception as e:
                                    print(f"[WARNING] Could not extract field '{key}': {str(e)}")
                                    continue
                    
                    col_meta_found = True
                    break
            
            if not col_meta_found:
                print("[ERROR] No column metadata found in any expected locations")
                print("[INFO] Available paths in file:")
                for path in sorted(all_paths):
                    print(f"  {path}")
                
                # Try to find any path containing "col" or "meta"
                candidate_paths = [p for p in all_paths if 'col' in p.lower() or 'meta' in p.lower()]
                if candidate_paths:
                    print(f"[INFO] Potential metadata paths found: {candidate_paths}")
                    print("[INFO] Try running in explore mode to see the structure")
                
                return False
            
            if not col_meta_data:
                print("[ERROR] No metadata fields could be extracted")
                print(f"[DEBUG] Metadata location: {successful_path}")
                print(f"[DEBUG] Available keys: {list(gctx_file[successful_path].keys()) if successful_path else 'None'}")
                return False
            
            # Create DataFrame from metadata
            try:
                # Ensure all arrays have the same length
                lengths = [len(v) for v in col_meta_data.values()]
                print(f"[DEBUG] Field lengths: {dict(zip(col_meta_data.keys(), lengths))}")
                
                if len(set(lengths)) > 1:
                    print(f"[WARNING] Metadata fields have different lengths!")
                    # Use the minimum length to avoid errors
                    min_length = min(lengths)
                    col_meta_data = {k: v[:min_length] for k, v in col_meta_data.items()}
                    print(f"[INFO] Truncated all fields to length {min_length}")
                
                # Create DataFrame
                col_meta_df = pd.DataFrame(col_meta_data)
                
                # Ensure we have the critical columns for drug analysis
                required_cols = ['pert_id', 'pert_iname']
                missing_cols = [col for col in required_cols if col not in col_meta_df.columns]
                if missing_cols:
                    print(f"[WARNING] Missing critical columns for drug analysis: {missing_cols}")
                    print(f"[INFO] Available columns: {list(col_meta_df.columns)}")
                
                # Write to TSV
                col_meta_df.to_csv(output_path, sep="\t", index=False)
                
                print(f"[SUCCESS] Saved col_meta to {output_path}")
                print(f"[INFO] Columns: {list(col_meta_df.columns)}")
                print(f"[INFO] Rows: {len(col_meta_df)}")
                
                # Show sample of the data
                print(f"\n=== Sample of extracted metadata ===")
                print(col_meta_df.head())
                
                # Show summary of key columns if they exist
                if 'pert_id' in col_meta_df.columns:
                    non_null_pert_ids = col_meta_df['pert_id'].notna().sum()
                    print(f"[INFO] Non-null pert_id entries: {non_null_pert_ids}/{len(col_meta_df)}")
                
                if 'pert_iname' in col_meta_df.columns:
                    non_null_names = col_meta_df['pert_iname'].notna().sum()
                    print(f"[INFO] Non-null pert_iname entries: {non_null_names}/{len(col_meta_df)}")
                
                return True
                
            except Exception as e:
                print(f"[ERROR] Failed to create DataFrame: {str(e)}")
                print(f"[DEBUG] Available data keys: {list(col_meta_data.keys())}")
                print(f"[DEBUG] Data types: {[(k, type(v[0]) if len(v) > 0 else 'empty') for k, v in col_meta_data.items()]}")
                return False
    
    except MemoryError:
        print(f"[ERROR] Not enough memory to read metadata. Try reducing max_memory_gb.")
        return False
    except Exception as e:
        print(f"[ERROR] Failed to read GCTX file: {str(e)}")
        import traceback
        traceback.print_exc()
        return False
# Add this function after the existing functions

def parse_signature_ids_to_metadata(input_file="lincs_col_meta.tsv", output_file="lincs_col_meta_parsed.tsv"):
    """
    Parse LINCS signature IDs to extract perturbagen IDs and other metadata
    """
    print("[INFO] Parsing signature IDs to extract metadata...")
    
    df = pd.read_csv(input_file, sep="\t")
    
    if 'id' not in df.columns:
        print("[ERROR] No 'id' column found in metadata file")
        return False
    
    # Parse signature IDs - format: {plate}_{cell}_{time}:{pert_id}-{batch}:{dose}
    parsed_data = []
    
    for idx, sig_id in enumerate(df['id']):
        try:
            # Split on colon first
            parts = str(sig_id).split(':')
            
            if len(parts) >= 3:
                # parts[0] = plate_cell_time
                # parts[1] = pert_id-batch_info  
                # parts[2] = dose
                
                plate_cell_time = parts[0]
                pert_batch = parts[1]
                dose = parts[2] if len(parts) > 2 else None
                
                # Extract pert_id (before the first dash in pert_batch)
                if '-' in pert_batch:
                    pert_id = pert_batch.split('-')[0]
                else:
                    pert_id = pert_batch
                
                # Extract plate, cell, time from first part
                pct_parts = plate_cell_time.split('_')
                if len(pct_parts) >= 3:
                    plate = pct_parts[0]
                    cell_id = pct_parts[1] 
                    pert_time = pct_parts[2]
                else:
                    plate = plate_cell_time
                    cell_id = None
                    pert_time = None
                
                parsed_data.append({
                    'id': sig_id,
                    'pert_id': pert_id,
                    'cell_id': cell_id,
                    'pert_time': pert_time,
                    'pert_dose': dose,
                    'plate': plate,
                    'full_pert_batch': pert_batch
                })
                
            else:
                # Couldn't parse - keep original
                parsed_data.append({
                    'id': sig_id,
                    'pert_id': None,
                    'cell_id': None, 
                    'pert_time': None,
                    'pert_dose': None,
                    'plate': None,
                    'full_pert_batch': None
                })
                
        except Exception as e:
            print(f"[WARNING] Could not parse signature ID '{sig_id}': {e}")
            parsed_data.append({
                'id': sig_id,
                'pert_id': None,
                'cell_id': None,
                'pert_time': None, 
                'pert_dose': None,
                'plate': None,
                'full_pert_batch': None
            })
    
    # Create new dataframe
    parsed_df = pd.DataFrame(parsed_data)
    
    # Save parsed metadata
    parsed_df.to_csv(output_file, sep="\t", index=False)
    
    print(f"[SUCCESS] Saved parsed metadata to {output_file}")
    print(f"[INFO] Columns: {list(parsed_df.columns)}")
    print(f"[INFO] Rows: {len(parsed_df)}")
    
    # Show summary
    print(f"\n=== Parsed Metadata Sample ===")
    print(parsed_df.head(10))
    
    # Show pert_id stats
    non_null_pert = parsed_df['pert_id'].notna().sum()
    unique_perts = parsed_df['pert_id'].nunique()
    print(f"\n[INFO] Successfully parsed pert_id from {non_null_pert}/{len(parsed_df)} signatures")
    print(f"[INFO] Found {unique_perts} unique perturbagens")
    
    return True

# Update the main function to include parsing option
def main():
    parser = argparse.ArgumentParser(description="Safely read GCTX files and extract metadata")
    parser.add_argument("--mode", choices=["explore", "extract", "parse"], default="explore",
                       help="Mode: 'explore' to view file structure, 'extract' to extract metadata, 'parse' to parse signature IDs")
    parser.add_argument("--gctx", default=None,
                       help="Path to GCTX file")
    parser.add_argument("--out", default="lincs_col_meta.tsv",
                       help="Output TSV file for metadata")
    parser.add_argument("--lines", type=int, default=10,
                       help="Number of lines to show in explore mode")
    parser.add_argument("--memory", type=float, default=2.0,
                       help="Maximum memory to use in GB")
    
    args = parser.parse_args()
    
    if args.mode == "explore":
        print(f"Exploring GCTX file: {args.gctx}")
        safe_read_gctx_top_lines(args.gctx, n_lines=args.lines, max_memory_gb=args.memory)
    
    elif args.mode == "extract":
        print(f"Extracting metadata from GCTX file: {args.gctx}")
        success = extract_lincs_col_metadata(args.gctx, args.out, max_memory_gb=args.memory)
        if success:
            print(f"[SUCCESS] Metadata extracted successfully!")
        else:
            print(f"[FAILED] Could not extract metadata")
            print(f"[TIP] Try running with --mode explore first to see the file structure")
    
    elif args.mode == "parse":
        print("Parsing signature IDs to extract perturbagen information...")
        success = parse_signature_ids_to_metadata(args.out, args.out.replace('.tsv', '_parsed.tsv'))
        if success:
            print(f"[SUCCESS] Signature IDs parsed successfully!")
            print(f"[INFO] Now use the '_parsed.tsv' file for drug analysis")
        else:
            print(f"[FAILED] Could not parse signature IDs")


if __name__ == "__main__":
    main()