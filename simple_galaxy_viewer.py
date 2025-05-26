#!/usr/bin/env python3
"""
Perfect Galaxy Combinations
Focus on perfecting the best filter combinations with fine-tuning options
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

def fine_tune_stretch(data, method='zscale', **kwargs):
    """Fine-tuned stretching with adjustable parameters"""
    if data is None or np.max(data) == 0:
        return np.zeros_like(data) if data is not None else None
    
    data_pos = np.maximum(data, 0)
    
    if method == 'zscale':
        low_cut = kwargs.get('low_cut', 1)
        high_cut = kwargs.get('high_cut', 99)
        vmin = np.percentile(data_pos, low_cut)
        vmax = np.percentile(data_pos, high_cut) 
        
        if vmax == vmin:
            return np.zeros_like(data_pos)
        
        stretched = (data_pos - vmin) / (vmax - vmin)
        stretched = np.clip(stretched, 0, 1)
        
    elif method == 'asinh_tuned':
        alpha = kwargs.get('alpha', 0.1)
        beta = kwargs.get('beta', 0.05)
        data_norm = data_pos / np.percentile(data_pos, 98)
        stretched = np.arcsinh(alpha * data_norm) / np.arcsinh(alpha) * (1 + beta)
        stretched = np.clip(stretched, 0, 1)
        
    elif method == 'gamma':
        gamma = kwargs.get('gamma', 0.6)
        data_norm = data_pos / np.percentile(data_pos, 99.5)
        stretched = np.power(np.clip(data_norm, 0, 1), gamma)
        
    else:  # linear
        stretched = data_pos / np.percentile(data_pos, 99.5)
        stretched = np.clip(stretched, 0, 1)
    
    return stretched

def create_perfect_rgb(r_data, g_data, b_data, 
                      r_params=None, g_params=None, b_params=None,
                      color_balance='auto', contrast=1.0, brightness=1.0, 
                      saturation=1.0, unsharp_strength=0.0):
    """Create perfectly tuned RGB composite"""
    
    # Default parameters
    if r_params is None: r_params = {'method': 'zscale', 'low_cut': 1, 'high_cut': 99}
    if g_params is None: g_params = {'method': 'zscale', 'low_cut': 1, 'high_cut': 99}
    if b_params is None: b_params = {'method': 'zscale', 'low_cut': 1, 'high_cut': 99}
    
    # Process each channel
    r_stretched = fine_tune_stretch(r_data, **r_params)
    g_stretched = fine_tune_stretch(g_data, **g_params)
    b_stretched = fine_tune_stretch(b_data, **b_params)
    
    if any(x is None for x in [r_stretched, g_stretched, b_stretched]):
        return None
    
    # Combine into RGB
    rgb = np.zeros((r_stretched.shape[0], r_stretched.shape[1], 3))
    rgb[:,:,0] = r_stretched
    rgb[:,:,1] = g_stretched
    rgb[:,:,2] = b_stretched
    
    # Color balance
    if color_balance == 'auto':
        # Balance based on image statistics
        for i in range(3):
            channel = rgb[:,:,i]
            if np.max(channel) > 0:
                # Normalize to similar 95th percentiles
                rgb[:,:,i] = channel / np.percentile(channel, 95)
    elif color_balance == 'manual':
        # Manual balance weights (could be parameters)
        weights = [1.0, 1.1, 1.2]  # Slightly boost green and blue
        for i in range(3):
            rgb[:,:,i] *= weights[i]
    
    # Apply unsharp masking if requested
    if unsharp_strength > 0:
        for i in range(3):
            channel = rgb[:,:,i]
            blurred = ndimage.gaussian_filter(channel, 1.5)
            rgb[:,:,i] = channel + unsharp_strength * (channel - blurred)
    
    # Global adjustments
    rgb = rgb * brightness * contrast
    
    # Saturation adjustment
    if saturation != 1.0:
        gray = 0.299 * rgb[:,:,0] + 0.587 * rgb[:,:,1] + 0.114 * rgb[:,:,2]
        for i in range(3):
            rgb[:,:,i] = gray + saturation * (rgb[:,:,i] - gray)
    
    return np.clip(rgb, 0, 1)

def create_masterpiece_versions(images, galaxy_name, combination_names):
    """Create multiple high-quality versions of the best combinations"""
    
    combinations = {
        'UV-V-I': ['f275w', 'f555w', 'f814w'],
        'B-V-I': ['f438w', 'f555w', 'f814w']
    }
    
    # Create figure for comparison
    fig, axes = plt.subplots(3, 2, figsize=(16, 24))
    fig.suptitle(f'{galaxy_name.upper()} - Masterpiece Versions', fontsize=20, fontweight='bold')
    
    saved_files = []
    
    for col, (combo_name, filters) in enumerate(combinations.items()):
        if not all(f in images for f in filters):
            continue
            
        print(f"Creating masterpiece versions for {combo_name}...")
        
        # Version 1: Clean and balanced
        rgb_clean = create_perfect_rgb(
            images[filters[0]], images[filters[1]], images[filters[2]],
            r_params={'method': 'zscale', 'low_cut': 0.5, 'high_cut': 99.5},
            g_params={'method': 'zscale', 'low_cut': 1, 'high_cut': 99},
            b_params={'method': 'zscale', 'low_cut': 1.5, 'high_cut': 98.5},
            color_balance='auto',
            brightness=1.0,
            saturation=1.1
        )
        
        axes[0, col].imshow(rgb_clean, origin='lower')
        axes[0, col].set_title(f'{combo_name} - Clean & Balanced', fontsize=14, fontweight='bold')
        axes[0, col].axis('off')
        
        # Version 2: Enhanced detail
        rgb_enhanced = create_perfect_rgb(
            images[filters[0]], images[filters[1]], images[filters[2]],
            r_params={'method': 'gamma', 'gamma': 0.6},
            g_params={'method': 'gamma', 'gamma': 0.65},
            b_params={'method': 'gamma', 'gamma': 0.7},
            color_balance='auto',
            brightness=1.1,
            saturation=1.3,
            unsharp_strength=0.2
        )
        
        axes[1, col].imshow(rgb_enhanced, origin='lower')
        axes[1, col].set_title(f'{combo_name} - Enhanced Detail', fontsize=14, fontweight='bold')
        axes[1, col].axis('off')
        
        # Version 3: Dramatic
        rgb_dramatic = create_perfect_rgb(
            images[filters[0]], images[filters[1]], images[filters[2]],
            r_params={'method': 'asinh_tuned', 'alpha': 0.2, 'beta': 0.1},
            g_params={'method': 'asinh_tuned', 'alpha': 0.15, 'beta': 0.05},
            b_params={'method': 'asinh_tuned', 'alpha': 0.1, 'beta': 0.0},
            color_balance='manual',
            brightness=1.2,
            saturation=1.5,
            unsharp_strength=0.3
        )
        
        axes[2, col].imshow(rgb_dramatic, origin='lower')
        axes[2, col].set_title(f'{combo_name} - Dramatic', fontsize=14, fontweight='bold')
        axes[2, col].axis('off')
        
        # Save high-resolution versions
        versions = [
            (rgb_clean, 'clean'),
            (rgb_enhanced, 'enhanced'), 
            (rgb_dramatic, 'dramatic')
        ]
        
        for rgb_data, version_name in versions:
            if rgb_data is not None:
                # Create individual high-res figure
                plt.figure(figsize=(15, 15))
                plt.imshow(rgb_data, origin='lower')
                plt.title(f'{galaxy_name.upper()} - {combo_name} ({version_name.title()})', 
                         fontsize=18, fontweight='bold', color='white', pad=20)
                plt.axis('off')
                plt.tight_layout()
                
                filename = f"{galaxy_name}_{combo_name.replace('-', '_')}_{version_name}.png"
                plt.savefig(filename, dpi=300, bbox_inches='tight', 
                           facecolor='black', edgecolor='none')
                plt.close()
                
                saved_files.append(filename)
                print(f"💾 Saved: {filename}")
    
    plt.tight_layout()
    plt.show()
    
    return fig, saved_files

def create_fine_tuned_gallery(images, galaxy_name):
    """Create gallery with fine-tuned parameters"""
    
    # Focus on the two best combinations
    combinations = [
        ('UV-V-I Classic', ['f275w', 'f555w', 'f814w'], 'Classic balanced look'),
        ('UV-V-I Enhanced', ['f275w', 'f555w', 'f814w'], 'Star formation emphasis'),
        ('B-V-I Classic', ['f438w', 'f555w', 'f814w'], 'Traditional optical'),
        ('B-V-I Enhanced', ['f438w', 'f555w', 'f814w'], 'Stellar population focus')
    ]
    
    fig, axes = plt.subplots(2, 2, figsize=(16, 16))
    fig.suptitle(f'{galaxy_name.upper()} - Fine-Tuned Gallery', fontsize=18, fontweight='bold')
    
    axes = axes.flatten()
    
    settings = [
        # UV-V-I Classic
        {'brightness': 1.0, 'saturation': 1.1, 'contrast': 1.0, 'unsharp_strength': 0.0},
        # UV-V-I Enhanced  
        {'brightness': 1.1, 'saturation': 1.4, 'contrast': 1.1, 'unsharp_strength': 0.2},
        # B-V-I Classic
        {'brightness': 1.0, 'saturation': 1.0, 'contrast': 1.0, 'unsharp_strength': 0.0},
        # B-V-I Enhanced
        {'brightness': 1.05, 'saturation': 1.2, 'contrast': 1.05, 'unsharp_strength': 0.15}
    ]
    
    for i, ((name, filters, desc), setting) in enumerate(zip(combinations, settings)):
        if not all(f in images for f in filters):
            continue
            
        rgb = create_perfect_rgb(
            images[filters[0]], images[filters[1]], images[filters[2]],
            color_balance='auto',
            **setting
        )
        
        if rgb is not None:
            axes[i].imshow(rgb, origin='lower')
            axes[i].set_title(f'{name}\n{desc}', fontsize=12, fontweight='bold')
            axes[i].axis('off')
    
    plt.tight_layout()
    plt.show()
    
    return fig

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
        print("Usage: python perfect_galaxy_combos.py /path/to/fits/directory")
        sys.exit(1)
    
    data_directory = sys.argv[1]
    
    print("🎨 Perfect Galaxy Combinations")
    print("="*50)
    print("Focusing on UV-V-I and B-V-I combinations")
    
    galaxy_data = find_phangs_files(data_directory)
    
    if not galaxy_data:
        print("No PHANGS data found!")
        sys.exit(1)
    
    for galaxy_name, galaxy_filters in galaxy_data.items():
        print(f"\n🌌 Perfecting {galaxy_name.upper()}")
        
        # Load images
        images = {}
        for filt in galaxy_filters:
            for file_type, file_path in galaxy_filters[filt].items():
                if 'drc' in file_type or 'sci' in file_type:
                    data = load_fits_image(file_path)
                    if data is not None:
                        images[filt] = data
                    break
        
        print(f"Loaded filters: {list(images.keys())}")
        
        # Create masterpiece versions
        print("Creating masterpiece versions...")
        fig1, saved_files = create_masterpiece_versions(images, galaxy_name, ['UV-V-I', 'B-V-I'])
        
        # Create fine-tuned gallery
        print("Creating fine-tuned gallery...")
        fig2 = create_fine_tuned_gallery(images, galaxy_name)
        
        print(f"\n✅ Created {len(saved_files)} high-resolution images:")
        for filename in saved_files:
            print(f"   💾 {filename}")
    
    print("\n🌟 Perfect combinations complete!")

if __name__ == "__main__":
    main()