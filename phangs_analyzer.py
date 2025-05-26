#!/usr/bin/env python3
"""
PHANGS Galaxy FITS File Analyzer
Basic tool to load, enhance, and visualize galaxy data from PHANGS survey

Usage: python phangs_analyzer.py /path/to/fits/directory
"""

import sys
import os
import glob
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import warnings
warnings.filterwarnings('ignore')

def find_phangs_files(directory):
    """
    Find and categorize PHANGS FITS files in directory
    Expected format: hlsp_phangs-hst_hst_wfc3-uvis_ngc4535_f275w_v1_err-drc-wht.fits
    """
    pattern = os.path.join(directory, "hlsp_phangs*.fits")
    files = glob.glob(pattern)
    
    if not files:
        print(f"No PHANGS FITS files found in {directory}")
        print("Expected format: hlsp_phangs-hst_hst_*_ngc*_f*w_v*_*.fits")
        return {}
    
    galaxy_data = {}
    
    print(f"Found {len(files)} FITS files")
    
    for file_path in files:
        filename = os.path.basename(file_path)
        print(f"Processing: {filename}")
        parts = filename.split('_')
        
        try:
            # Parse filename parts for your specific format
            # hlsp_phangs-hst_hst_wfc3-uvis_ngc4535_f275w_v1_err-drc-wht.fits
            #   0      1         2      3        4       5  6      7
            
            galaxy = parts[4]  # e.g., ngc4535
            filter_name = parts[5]  # e.g., f275w
            file_type = parts[-1].replace('.fits', '')  # e.g., err-drc-wht
            
            print(f"  -> Galaxy: {galaxy}, Filter: {filter_name}, Type: {file_type}")
            
            if galaxy not in galaxy_data:
                galaxy_data[galaxy] = {}
            if filter_name not in galaxy_data[galaxy]:
                galaxy_data[galaxy][filter_name] = {}
            
            galaxy_data[galaxy][filter_name][file_type] = file_path
            
        except IndexError:
            print(f"Skipping file with unexpected format: {filename}")
            print(f"  Parts: {parts}")
    
    return galaxy_data

def load_fits_image(file_path):
    """Load FITS image data, handling common formats"""
    try:
        with fits.open(file_path) as hdul:
            # Try different extensions
            if len(hdul) > 1 and hdul[1].data is not None:
                data = hdul[1].data
                header = hdul[1].header
            else:
                data = hdul[0].data
                header = hdul[0].header
            
            # Clean up data
            data = np.nan_to_num(data, nan=0.0, posinf=0.0, neginf=0.0)
            
            return data, header
            
    except Exception as e:
        print(f"Error loading {file_path}: {e}")
        return None, None

def enhance_image(data, method='asinh', clip_percent=99.5):
    """
    Basic image enhancement for astronomy data
    
    Parameters:
    - data: 2D numpy array
    - method: 'asinh', 'sqrt', 'log'
    - clip_percent: percentile for clipping bright pixels
    """
    if data is None:
        return None
    
    # Remove negative values and clip extreme values
    data_pos = np.maximum(data, 0)
    clip_val = np.percentile(data_pos[data_pos > 0], clip_percent)
    data_clipped = np.minimum(data_pos, clip_val)
    
    if method == 'asinh':
        # Asinh stretch - good for wide dynamic range
        enhanced = np.arcsinh(data_clipped / (clip_val * 0.1))
    elif method == 'sqrt':
        # Square root stretch
        enhanced = np.sqrt(data_clipped / clip_val)
    elif method == 'log':
        # Log stretch
        enhanced = np.log10(data_clipped / clip_val + 1)
    else:
        enhanced = data_clipped / clip_val
    
    return enhanced

def plot_galaxy_overview(galaxy_name, galaxy_filters, max_filters=6):
    """Create overview plot showing all available filters for a galaxy"""
    
    # Get science images (files containing 'drc' or 'sci')
    science_images = {}
    for filt, file_types in galaxy_filters.items():
        science_file = None
        
        # Look for science data files (prefer drc, then sci)
        for file_type, file_path in file_types.items():
            if 'drc' in file_type or 'sci' in file_type:
                science_file = file_path
                break
        
        if science_file:
            data, header = load_fits_image(science_file)
            if data is not None:
                science_images[filt] = data
                print(f"  Loaded {filt}: {data.shape}, range: {data.min():.2e} to {data.max():.2e}")
    
    if not science_images:
        print(f"No science images found for {galaxy_name}")
        print(f"Available file types: {[list(f.keys()) for f in galaxy_filters.values()]}")
        return
    
    # Limit number of filters to display
    filters_to_plot = list(science_images.keys())[:max_filters]
    n_filters = len(filters_to_plot)
    
    # Set up plot grid
    if n_filters <= 3:
        rows, cols = 1, n_filters
        figsize = (5 * n_filters, 5)
    else:
        rows = 2
        cols = (n_filters + 1) // 2
        figsize = (5 * cols, 10)
    
    fig, axes = plt.subplots(rows, cols, figsize=figsize)
    if n_filters == 1:
        axes = [axes]
    elif rows == 1:
        axes = axes
    else:
        axes = axes.flatten()
    
    fig.suptitle(f'{galaxy_name.upper()} - Multi-Filter Overview', fontsize=16, fontweight='bold')
    
    for i, filt in enumerate(filters_to_plot):
        data = science_images[filt]
        enhanced = enhance_image(data, 'asinh')
        
        if enhanced is not None:
            im = axes[i].imshow(enhanced, cmap='viridis', origin='lower')
            axes[i].set_title(f'{filt.upper()}', fontsize=12, fontweight='bold')
            axes[i].axis('off')
            
            # Add colorbar
            plt.colorbar(im, ax=axes[i], shrink=0.8)
            
            # Add image info
            axes[i].text(0.02, 0.98, f'{data.shape[0]}×{data.shape[1]}', 
                        transform=axes[i].transAxes, fontsize=8, color='white',
                        verticalalignment='top', bbox=dict(boxstyle='round', facecolor='black', alpha=0.5))
    
    # Hide unused subplots
    for i in range(n_filters, len(axes)):
        axes[i].axis('off')
    
    plt.tight_layout()
    plt.show()
    
    return fig

def create_enhanced_comparison(galaxy_name, galaxy_filters, filter_name=None):
    """Show different enhancement methods for a single filter"""
    
    # Choose filter to analyze
    if filter_name is None:
        available_filters = []
        for filt, file_types in galaxy_filters.items():
            for file_type in file_types.keys():
                if 'drc' in file_type or 'sci' in file_type:
                    available_filters.append(filt)
                    break
        
        if not available_filters:
            print("No science images available")
            print(f"Available file types: {[list(f.keys()) for f in galaxy_filters.values()]}")
            return
        filter_name = available_filters[0]  # Use first available
    
    # Find science file for this filter
    science_file = None
    if filter_name in galaxy_filters:
        for file_type, file_path in galaxy_filters[filter_name].items():
            if 'drc' in file_type or 'sci' in file_type:
                science_file = file_path
                break
    
    if science_file is None:
        print(f"No science file found for filter {filter_name}")
        return
    
    # Load data
    data, header = load_fits_image(science_file)
    if data is None:
        return
    
    print(f"Loaded {filter_name}: {data.shape}, range: {data.min():.2e} to {data.max():.2e}")
    
    # Create different enhancements
    methods = ['asinh', 'sqrt', 'log']
    enhanced_images = {}
    
    for method in methods:
        enhanced_images[method] = enhance_image(data, method)
    
    # Plot comparison
    fig, axes = plt.subplots(1, 4, figsize=(20, 5))
    fig.suptitle(f'{galaxy_name.upper()} - {filter_name.upper()} Enhancement Comparison', 
                 fontsize=16, fontweight='bold')
    
    # Original
    im0 = axes[0].imshow(data, cmap='gray', origin='lower')
    axes[0].set_title('Original')
    axes[0].axis('off')
    plt.colorbar(im0, ax=axes[0], shrink=0.8)
    
    # Enhanced versions
    cmaps = ['viridis', 'plasma', 'inferno']
    for i, (method, cmap) in enumerate(zip(methods, cmaps)):
        enhanced = enhanced_images[method]
        if enhanced is not None:
            im = axes[i+1].imshow(enhanced, cmap=cmap, origin='lower')
            axes[i+1].set_title(f'{method.capitalize()} Stretch')
            axes[i+1].axis('off')
            plt.colorbar(im, ax=axes[i+1], shrink=0.8)
    
    plt.tight_layout()
    plt.show()
    
    return fig

def print_data_summary(galaxy_data):
    """Print summary of available data"""
    print("\n" + "="*60)
    print("PHANGS DATA SUMMARY")
    print("="*60)
    
    for galaxy, filters in galaxy_data.items():
        print(f"\n📍 {galaxy.upper()}")
        print(f"   Available filters: {len(filters)}")
        
        for filt, file_types in filters.items():
            print(f"   • {filt.upper()}: {', '.join(file_types.keys())}")
            
            # Show image dimensions for science files
            science_file = None
            for file_type, file_path in file_types.items():
                if 'drc' in file_type or 'sci' in file_type:
                    science_file = file_path
                    break
            
            if science_file:
                data, _ = load_fits_image(science_file)
                if data is not None:
                    print(f"     └─ Size: {data.shape[0]} × {data.shape[1]} pixels")

def main():
    """Main analysis function"""
    if len(sys.argv) != 2:
        print("Usage: python phangs_analyzer.py /path/to/fits/directory")
        sys.exit(1)
    
    data_directory = sys.argv[1]
    
    if not os.path.exists(data_directory):
        print(f"Directory not found: {data_directory}")
        sys.exit(1)
    
    print(f"🔍 Searching for PHANGS FITS files in: {data_directory}")
    
    # Find and categorize files
    galaxy_data = find_phangs_files(data_directory)
    
    if not galaxy_data:
        print("No PHANGS data found!")
        sys.exit(1)
    
    # Print summary
    print_data_summary(galaxy_data)
    
    # Analyze each galaxy
    for galaxy_name, galaxy_filters in galaxy_data.items():
        print(f"\n🌌 Analyzing {galaxy_name.upper()}...")
        
        # Create overview plot
        print("   Creating multi-filter overview...")
        plot_galaxy_overview(galaxy_name, galaxy_filters)
        
        # Create enhancement comparison
        print("   Creating enhancement comparison...")
        create_enhanced_comparison(galaxy_name, galaxy_filters)
    
    print("\n✅ Analysis complete!")
    print("\nNext steps:")
    print("• Try different filters in create_enhanced_comparison()")
    print("• Experiment with enhancement parameters") 
    print("• Add your own analysis functions")

if __name__ == "__main__":
    main()