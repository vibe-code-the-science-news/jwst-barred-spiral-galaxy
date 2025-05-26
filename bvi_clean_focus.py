#!/usr/bin/env python3
"""
Focused B-V-I Image Creator
Perfect the B-V-I combination that looked great, avoiding display issues
"""

import sys
import os
import glob
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy import ndimage
import warnings
warnings.filterwarnings('ignore')

def load_fits_image(file_path):
    """Load FITS image data"""
    try:
        with fits.open(file_path) as hdul:
            if len(hdul) > 1 and hdul[1].data is not None:
                data = hdul[1].data
            else:
                data = hdul[0].data
            
            data = np.nan_to_num(data, nan=0.0, posinf=0.0, neginf=0.0)
            return data
    except Exception as e:
        print(f"Error loading {file_path}: {e}")
        return None

def perfect_stretch(data, method='balanced_zscale', **kwargs):
    """Perfect stretching for the B-V-I combination"""
    if data is None or np.max(data) == 0:
        return np.zeros_like(data) if data is not None else None
    
    data_pos = np.maximum(data, 0)
    
    if method == 'balanced_zscale':
        # Optimized for clean, balanced look
        low_cut = kwargs.get('low_cut', 0.5)
        high_cut = kwargs.get('high_cut', 99.2)
        vmin = np.percentile(data_pos, low_cut)
        vmax = np.percentile(data_pos, high_cut)
        
        if vmax == vmin:
            return np.zeros_like(data_pos)
        
        stretched = (data_pos - vmin) / (vmax - vmin)
        stretched = np.clip(stretched, 0, 1)
        
    elif method == 'smooth_gamma':
        # Smooth gamma correction
        gamma = kwargs.get('gamma', 0.65)
        data_norm = data_pos / np.percentile(data_pos, 99.5)
        stretched = np.power(np.clip(data_norm, 0, 1), gamma)
        
    elif method == 'gentle_asinh':
        # Very gentle asinh for subtle enhancement
        alpha = kwargs.get('alpha', 0.05)
        data_norm = data_pos / np.percentile(data_pos, 98)
        stretched = np.arcsinh(alpha * data_norm) / np.arcsinh(alpha)
        stretched = np.clip(stretched, 0, 1)
        
    else:  # linear fallback
        stretched = data_pos / np.percentile(data_pos, 99.5)
        stretched = np.clip(stretched, 0, 1)
    
    return stretched

def create_perfect_bvi(b_data, v_data, i_data, style='clean'):
    """Create perfect B-V-I composite"""
    
    if style == 'clean':
        # Clean and balanced - what looked great!
        b_stretched = perfect_stretch(b_data, 'balanced_zscale', low_cut=0.3, high_cut=99.5)
        v_stretched = perfect_stretch(v_data, 'balanced_zscale', low_cut=0.5, high_cut=99.2)
        i_stretched = perfect_stretch(i_data, 'balanced_zscale', low_cut=0.8, high_cut=98.8)
        
        brightness = 1.0
        saturation = 1.05
        unsharp = 0.0
        
    elif style == 'enhanced':
        # Gentle enhancement
        b_stretched = perfect_stretch(b_data, 'smooth_gamma', gamma=0.6)
        v_stretched = perfect_stretch(v_data, 'smooth_gamma', gamma=0.65)
        i_stretched = perfect_stretch(i_data, 'smooth_gamma', gamma=0.7)
        
        brightness = 1.05
        saturation = 1.15
        unsharp = 0.1
        
    elif style == 'detailed':
        # More detail oriented
        b_stretched = perfect_stretch(b_data, 'gentle_asinh', alpha=0.08)
        v_stretched = perfect_stretch(v_data, 'gentle_asinh', alpha=0.06)
        i_stretched = perfect_stretch(i_data, 'gentle_asinh', alpha=0.04)
        
        brightness = 1.1
        saturation = 1.25
        unsharp = 0.2
        
    else:
        return None
    
    if any(x is None for x in [b_stretched, v_stretched, i_stretched]):
        return None
    
    # Combine into RGB
    rgb = np.zeros((b_stretched.shape[0], b_stretched.shape[1], 3))
    rgb[:,:,0] = b_stretched  # Red channel
    rgb[:,:,1] = v_stretched  # Green channel  
    rgb[:,:,2] = i_stretched  # Blue channel
    
    # Color balance - optimize for B-V-I
    for i in range(3):
        channel = rgb[:,:,i]
        if np.max(channel) > 0:
            # Use 96th percentile for B-V-I balance
            rgb[:,:,i] = channel / np.percentile(channel, 96)
    
    # Apply unsharp masking if requested
    if unsharp > 0:
        for i in range(3):
            channel = rgb[:,:,i]
            blurred = ndimage.gaussian_filter(channel, 1.2)
            rgb[:,:,i] = channel + unsharp * (channel - blurred)
    
    # Global adjustments
    rgb = rgb * brightness
    
    # Saturation adjustment optimized for stellar populations
    if saturation != 1.0:
        gray = 0.299 * rgb[:,:,0] + 0.587 * rgb[:,:,1] + 0.114 * rgb[:,:,2]
        for i in range(3):
            rgb[:,:,i] = gray + saturation * (rgb[:,:,i] - gray)
    
    return np.clip(rgb, 0, 1)

def create_bvi_masterpieces(images, galaxy_name):
    """Create B-V-I masterpieces without display issues"""
    
    if not all(f in images for f in ['f438w', 'f555w', 'f814w']):
        print("Missing required filters for B-V-I")
        return []
    
    styles = [
        ('clean', 'Clean & Balanced'),
        ('enhanced', 'Gently Enhanced'), 
        ('detailed', 'Detail Enhanced')
    ]
    
    saved_files = []
    
    print("Creating B-V-I masterpieces...")
    
    for style_code, style_name in styles:
        print(f"  → Creating {style_name} version...")
        
        rgb = create_perfect_bvi(
            images['f438w'],  # B-band
            images['f555w'],  # V-band  
            images['f814w'],  # I-band
            style=style_code
        )
        
        if rgb is not None:
            # Create high-resolution figure
            plt.figure(figsize=(16, 16), dpi=100)
            plt.imshow(rgb, origin='lower')
            plt.title(f'{galaxy_name.upper()} B-V-I Composite\n{style_name}', 
                     fontsize=20, fontweight='bold', color='white', pad=30)
            plt.axis('off')
            plt.tight_layout()
            
            # Save without showing (avoids display issues)
            filename = f"{galaxy_name}_BVI_{style_code}_perfected.png"
            plt.savefig(filename, dpi=300, bbox_inches='tight', 
                       facecolor='black', edgecolor='none', pad_inches=0.1)
            plt.close()  # Important: close figure to free memory
            
            saved_files.append(filename)
            print(f"    💾 Saved: {filename}")
    
    return saved_files

def create_bvi_comparison(images, galaxy_name):
    """Create comparison of B-V-I processing approaches"""
    
    if not all(f in images for f in ['f438w', 'f555w', 'f814w']):
        return None
    
    print("Creating B-V-I comparison...")
    
    # Create different versions
    versions = {
        'Original Clean': create_perfect_bvi(images['f438w'], images['f555w'], images['f814w'], 'clean'),
        'Optimized': create_perfect_bvi(images['f438w'], images['f555w'], images['f814w'], 'enhanced'),
        'Detail Focus': create_perfect_bvi(images['f438w'], images['f555w'], images['f814w'], 'detailed')
    }
    
    # Create comparison figure
    plt.figure(figsize=(18, 6))
    
    for i, (name, rgb) in enumerate(versions.items()):
        if rgb is not None:
            plt.subplot(1, 3, i + 1)
            plt.imshow(rgb, origin='lower')
            plt.title(f'{name}', fontsize=14, fontweight='bold')
            plt.axis('off')
    
    plt.suptitle(f'{galaxy_name.upper()} - B-V-I Processing Comparison', 
                 fontsize=16, fontweight='bold')
    plt.tight_layout()
    
    # Save comparison
    comparison_file = f"{galaxy_name}_BVI_comparison.png"
    plt.savefig(comparison_file, dpi=200, bbox_inches='tight', facecolor='white')
    plt.close()
    
    print(f"💾 Saved comparison: {comparison_file}")
    return comparison_file

def find_phangs_files(directory):
    """Find PHANGS files"""
    pattern = os.path.join(directory, "hlsp_phangs*.fits")
    files = glob.glob(pattern)
    
    galaxy_data = {}
    for file_path in files:
        filename = os.path.basename(file_path)
        parts = filename.split('_')
        
        try:
            galaxy = parts[4]
            filter_name = parts[5]
            file_type = parts[-1].replace('.fits', '')
            
            if galaxy not in galaxy_data:
                galaxy_data[galaxy] = {}
            if filter_name not in galaxy_data[galaxy]:
                galaxy_data[galaxy][filter_name] = {}
            
            galaxy_data[galaxy][filter_name][file_type] = file_path
        except IndexError:
            continue
    
    return galaxy_data

def main():
    """Main function"""
    if len(sys.argv) != 2:
        print("Usage: python focus_bvi_creator.py /path/to/fits/directory")
        sys.exit(1)
    
    data_directory = sys.argv[1]
    
    print("🌌 Focused B-V-I Image Creator")
    print("="*50)
    print("Perfecting the B-V-I combination that looked great!")
    
    galaxy_data = find_phangs_files(data_directory)
    
    if not galaxy_data:
        print("No PHANGS data found!")
        sys.exit(1)
    
    for galaxy_name, galaxy_filters in galaxy_data.items():
        print(f"\n🎨 Perfecting B-V-I for {galaxy_name.upper()}")
        
        # Load images
        images = {}
        for filt in galaxy_filters:
            for file_type, file_path in galaxy_filters[filt].items():
                if 'drc' in file_type or 'sci' in file_type:
                    data = load_fits_image(file_path)
                    if data is not None:
                        images[filt] = data
                        print(f"  Loaded {filt}: {data.shape}")
                    break
        
        # Create B-V-I masterpieces
        saved_files = create_bvi_masterpieces(images, galaxy_name)
        
        # Create comparison
        comparison_file = create_bvi_comparison(images, galaxy_name)
        if comparison_file:
            saved_files.append(comparison_file)
        
        print(f"\n✅ Created {len(saved_files)} B-V-I images:")
        for filename in saved_files:
            print(f"   🌟 {filename}")
    
    print(f"\n🎉 B-V-I perfection complete!")
    print("Focus on the 'clean' version - that's the one that looked great!")

if __name__ == "__main__":
    main()