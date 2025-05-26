#!/usr/bin/env python3
"""
Perfect UV-V-I Creator
Focus on perfecting the UV-V-I combination that looked amazing
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

def perfect_uv_stretch(data, method='uv_optimized', **kwargs):
    """Perfect stretching optimized for UV-V-I combination"""
    if data is None or np.max(data) == 0:
        return np.zeros_like(data) if data is not None else None
    
    data_pos = np.maximum(data, 0)
    
    if method == 'uv_optimized':
        # Optimized specifically for UV data (can be noisier)
        low_cut = kwargs.get('low_cut', 0.3)
        high_cut = kwargs.get('high_cut', 99.0)
        vmin = np.percentile(data_pos, low_cut)
        vmax = np.percentile(data_pos, high_cut)
        
        if vmax == vmin:
            return np.zeros_like(data_pos)
        
        stretched = (data_pos - vmin) / (vmax - vmin)
        stretched = np.clip(stretched, 0, 1)
        
    elif method == 'optical_balanced':
        # Optimized for V and I bands
        low_cut = kwargs.get('low_cut', 0.5)
        high_cut = kwargs.get('high_cut', 99.3)
        vmin = np.percentile(data_pos, low_cut)
        vmax = np.percentile(data_pos, high_cut)
        
        if vmax == vmin:
            return np.zeros_like(data_pos)
        
        stretched = (data_pos - vmin) / (vmax - vmin)
        stretched = np.clip(stretched, 0, 1)
        
    elif method == 'gentle_gamma':
        # Gentle gamma for enhanced version
        gamma = kwargs.get('gamma', 0.65)
        data_norm = data_pos / np.percentile(data_pos, 99.2)
        stretched = np.power(np.clip(data_norm, 0, 1), gamma)
        
    elif method == 'star_formation_focus':
        # Emphasize star formation regions (for UV)
        alpha = kwargs.get('alpha', 0.06)
        data_norm = data_pos / np.percentile(data_pos, 97)
        stretched = np.arcsinh(alpha * data_norm) / np.arcsinh(alpha)
        stretched = np.clip(stretched, 0, 1)
        
    else:  # linear fallback
        stretched = data_pos / np.percentile(data_pos, 99.5)
        stretched = np.clip(stretched, 0, 1)
    
    return stretched

def create_perfect_uv_v_i(uv_data, v_data, i_data, style='clean'):
    """Create perfect UV-V-I composite"""
    
    if style == 'clean':
        # Clean and balanced - the style you loved!
        uv_stretched = perfect_uv_stretch(uv_data, 'uv_optimized', low_cut=0.2, high_cut=99.0)
        v_stretched = perfect_uv_stretch(v_data, 'optical_balanced', low_cut=0.4, high_cut=99.4)
        i_stretched = perfect_uv_stretch(i_data, 'optical_balanced', low_cut=0.6, high_cut=99.1)
        
        brightness = 1.0
        saturation = 1.1
        unsharp = 0.0
        
    elif style == 'star_formation':
        # Emphasize star formation regions
        uv_stretched = perfect_uv_stretch(uv_data, 'star_formation_focus', alpha=0.08)
        v_stretched = perfect_uv_stretch(v_data, 'gentle_gamma', gamma=0.7)
        i_stretched = perfect_uv_stretch(i_data, 'gentle_gamma', gamma=0.75)
        
        brightness = 1.05
        saturation = 1.3
        unsharp = 0.1
        
    elif style == 'enhanced_clean':
        # Enhanced but still clean
        uv_stretched = perfect_uv_stretch(uv_data, 'uv_optimized', low_cut=0.1, high_cut=98.8)
        v_stretched = perfect_uv_stretch(v_data, 'gentle_gamma', gamma=0.68)
        i_stretched = perfect_uv_stretch(i_data, 'gentle_gamma', gamma=0.72)
        
        brightness = 1.03
        saturation = 1.2
        unsharp = 0.05
        
    elif style == 'structure_focus':
        # Focus on spiral structure
        uv_stretched = perfect_uv_stretch(uv_data, 'uv_optimized', low_cut=0.15, high_cut=98.5)
        v_stretched = perfect_uv_stretch(v_data, 'optical_balanced', low_cut=0.3, high_cut=99.2)
        i_stretched = perfect_uv_stretch(i_data, 'optical_balanced', low_cut=0.5, high_cut=98.9)
        
        brightness = 1.08
        saturation = 1.15
        unsharp = 0.15
        
    else:
        return None
    
    if any(x is None for x in [uv_stretched, v_stretched, i_stretched]):
        return None
    
    # Combine into RGB - UV-V-I mapping
    rgb = np.zeros((uv_stretched.shape[0], uv_stretched.shape[1], 3))
    rgb[:,:,0] = uv_stretched  # Red channel - UV (star formation)
    rgb[:,:,1] = v_stretched   # Green channel - V (overall structure)
    rgb[:,:,2] = i_stretched   # Blue channel - I (older stars)
    
    # Color balance optimized for UV-V-I
    # UV can be dimmer, so adjust balance carefully
    channel_weights = [0.95, 0.97, 0.98]  # Slightly boost UV relative to others
    
    for i in range(3):
        channel = rgb[:,:,i]
        if np.max(channel) > 0:
            # Use different percentiles for optimal UV-V-I balance
            percentile = [94, 96, 97][i]  # Lower for UV to boost it
            rgb[:,:,i] = (channel / np.percentile(channel, percentile)) * channel_weights[i]
    
    # Apply unsharp masking if requested
    if unsharp > 0:
        for i in range(3):
            channel = rgb[:,:,i]
            blurred = ndimage.gaussian_filter(channel, 1.0)
            rgb[:,:,i] = channel + unsharp * (channel - blurred)
    
    # Global adjustments
    rgb = rgb * brightness
    
    # Saturation adjustment optimized for star formation regions
    if saturation != 1.0:
        # Use luminance weights optimized for UV-V-I
        gray = 0.25 * rgb[:,:,0] + 0.65 * rgb[:,:,1] + 0.1 * rgb[:,:,2]
        for i in range(3):
            rgb[:,:,i] = gray + saturation * (rgb[:,:,i] - gray)
    
    return np.clip(rgb, 0, 1)

def create_uv_v_i_masterpieces(images, galaxy_name):
    """Create UV-V-I masterpieces"""
    
    if not all(f in images for f in ['f275w', 'f555w', 'f814w']):
        print("Missing required filters for UV-V-I")
        print(f"Available: {list(images.keys())}")
        print("Need: f275w (UV), f555w (V), f814w (I)")
        return []
    
    styles = [
        ('clean', 'Clean & Balanced - Your Favorite!'),
        ('enhanced_clean', 'Enhanced Clean'),
        ('star_formation', 'Star Formation Focus'),
        ('structure_focus', 'Spiral Structure Focus')
    ]
    
    saved_files = []
    
    print("Creating UV-V-I masterpieces...")
    print("UV (F275W) → Red channel (star formation regions)")
    print("V (F555W)  → Green channel (overall structure)")  
    print("I (F814W)  → Blue channel (older stellar populations)")
    
    for style_code, style_name in styles:
        print(f"  → Creating {style_name}...")
        
        rgb = create_perfect_uv_v_i(
            images['f275w'],  # UV-band → Red
            images['f555w'],  # V-band → Green
            images['f814w'],  # I-band → Blue
            style=style_code
        )
        
        if rgb is not None:
            # Create high-resolution figure
            plt.figure(figsize=(16, 16), dpi=100)
            plt.imshow(rgb, origin='lower')
            
            title_text = f'{galaxy_name.upper()} UV-V-I Composite\n{style_name}'
            if style_code == 'clean':
                title_text += '\n⭐ YOUR FAVORITE STYLE ⭐'
                
            plt.title(title_text, fontsize=20, fontweight='bold', 
                     color='white', pad=30)
            plt.axis('off')
            plt.tight_layout()
            
            # Save without showing
            filename = f"{galaxy_name}_UV_V_I_{style_code}_perfected.png"
            plt.savefig(filename, dpi=300, bbox_inches='tight', 
                       facecolor='black', edgecolor='none', pad_inches=0.1)
            plt.close()
            
            saved_files.append(filename)
            print(f"    💾 Saved: {filename}")
            
            # Extra special treatment for the clean version
            if style_code == 'clean':
                # Create an even larger version
                plt.figure(figsize=(20, 20), dpi=100)
                plt.imshow(rgb, origin='lower')
                plt.title(f'{galaxy_name.upper()} UV-V-I\nYour Perfect Clean Style', 
                         fontsize=24, fontweight='bold', color='white', pad=40)
                plt.axis('off')
                plt.tight_layout()
                
                large_filename = f"{galaxy_name}_UV_V_I_clean_LARGE.png"
                plt.savefig(large_filename, dpi=300, bbox_inches='tight', 
                           facecolor='black', edgecolor='none', pad_inches=0.15)
                plt.close()
                
                saved_files.append(large_filename)
                print(f"    🌟 Saved LARGE version: {large_filename}")
    
    return saved_files

def create_uv_v_i_comparison(images, galaxy_name):
    """Create UV-V-I comparison chart"""
    
    if not all(f in images for f in ['f275w', 'f555w', 'f814w']):
        return None
    
    print("Creating UV-V-I style comparison...")
    
    # Create different versions
    versions = {
        'Your Favorite Clean': create_perfect_uv_v_i(images['f275w'], images['f555w'], images['f814w'], 'clean'),
        'Enhanced Clean': create_perfect_uv_v_i(images['f275w'], images['f555w'], images['f814w'], 'enhanced_clean'),
        'Star Formation Focus': create_perfect_uv_v_i(images['f275w'], images['f555w'], images['f814w'], 'star_formation'),
        'Structure Focus': create_perfect_uv_v_i(images['f275w'], images['f555w'], images['f814w'], 'structure_focus')
    }
    
    # Create comparison figure
    plt.figure(figsize=(20, 10))
    
    for i, (name, rgb) in enumerate(versions.items()):
        if rgb is not None:
            plt.subplot(2, 2, i + 1)
            plt.imshow(rgb, origin='lower')
            title = f'{name}'
            if i == 0:  # First one is the favorite
                title += '\n⭐ YOUR CHOICE ⭐'
            plt.title(title, fontsize=14, fontweight='bold')
            plt.axis('off')
    
    plt.suptitle(f'{galaxy_name.upper()} - UV-V-I Style Comparison\nUV (Red) + V (Green) + I (Blue)', 
                 fontsize=18, fontweight='bold')
    plt.tight_layout()
    
    # Save comparison
    comparison_file = f"{galaxy_name}_UV_V_I_comparison.png"
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
        print("Usage: python perfect_uv_v_i_creator.py /path/to/fits/directory")
        sys.exit(1)
    
    data_directory = sys.argv[1]
    
    print("🌌 Perfect UV-V-I Creator")
    print("="*50)
    print("Perfecting the UV-V-I combination you love!")
    print("UV = Star formation regions (hot young stars)")
    print("V  = Overall galaxy structure") 
    print("I  = Older stellar populations")
    
    galaxy_data = find_phangs_files(data_directory)
    
    if not galaxy_data:
        print("No PHANGS data found!")
        sys.exit(1)
    
    for galaxy_name, galaxy_filters in galaxy_data.items():
        print(f"\n🎨 Perfecting UV-V-I for {galaxy_name.upper()}")
        
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
        
        # Create UV-V-I masterpieces
        saved_files = create_uv_v_i_masterpieces(images, galaxy_name)
        
        # Create comparison
        comparison_file = create_uv_v_i_comparison(images, galaxy_name)
        if comparison_file:
            saved_files.append(comparison_file)
        
        print(f"\n✅ Created {len(saved_files)} UV-V-I images:")
        for filename in saved_files:
            if 'clean' in filename.lower():
                print(f"   ⭐ {filename} ← YOUR FAVORITE!")
            else:
                print(f"   🌟 {filename}")
    
    print(f"\n🎉 UV-V-I perfection complete!")
    print("The 'clean' versions are optimized for your taste!")
    print("The LARGE version is extra high-resolution for printing/presentations!")

if __name__ == "__main__":
    main()