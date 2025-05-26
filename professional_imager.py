#!/usr/bin/env python3
"""
Professional Galaxy Image Creator
Creates publication-quality images like those seen in astronomy news articles

Usage: python professional_imager.py /path/to/fits/directory
"""

import sys
import os
import glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from astropy.io import fits
from scipy import ndimage, signal
from skimage import exposure, filters, restoration
import warnings
warnings.filterwarnings('ignore')

def load_fits_image(file_path):
    """Load FITS image data, handling common formats"""
    try:
        with fits.open(file_path) as hdul:
            if len(hdul) > 1 and hdul[1].data is not None:
                data = hdul[1].data
                header = hdul[1].header
            else:
                data = hdul[0].data
                header = hdul[0].header
            
            data = np.nan_to_num(data, nan=0.0, posinf=0.0, neginf=0.0)
            return data, header
    except Exception as e:
        print(f"Error loading {file_path}: {e}")
        return None, None

def professional_stretch(data, method='lupton', **kwargs):
    """
    Professional image stretching methods used in publication images
    """
    if data is None:
        return None
    
    # Clean data
    data_clean = np.maximum(data, 0)
    
    if method == 'lupton':
        # Lupton et al. (2004) asinh stretch - used by SDSS, HST
        alpha = kwargs.get('alpha', 1e-10)
        Q = kwargs.get('Q', 8)
        stretch = np.arcsinh(alpha * Q * data_clean) / Q
        
    elif method == 'zscale':
        # IRAF ZScale algorithm - good for astronomical images
        # Simplified version
        sorted_data = np.sort(data_clean.flatten())
        n = len(sorted_data)
        z1 = sorted_data[int(0.25 * n)]
        z2 = sorted_data[int(0.75 * n)]
        stretch = np.clip((data_clean - z1) / (z2 - z1), 0, 1)
        
    elif method == 'power':
        # Power law stretch
        gamma = kwargs.get('gamma', 0.5)
        stretch = np.power(data_clean / np.percentile(data_clean, 99.5), gamma)
        
    elif method == 'adaptive':
        # Adaptive histogram equalization
        # Normalize data to 0-1 range first
        data_norm = data_clean / np.percentile(data_clean, 99.9)
        data_norm = np.clip(data_norm, 0, 1)  # Ensure 0-1 range
        stretch = exposure.equalize_adapthist(data_norm, clip_limit=0.01)
        
    else:  # asinh fallback
        stretch = np.arcsinh(data_clean / np.percentile(data_clean, 95))
    
    return np.clip(stretch, 0, 1)

def unsharp_mask_pro(image, radius=2.0, amount=2.0, threshold=0.0):
    """
    Professional unsharp masking to enhance fine detail
    """
    # Create Gaussian blur
    blurred = ndimage.gaussian_filter(image, radius)
    
    # Create mask
    mask = image - blurred
    
    # Apply threshold
    mask = np.where(np.abs(mask) < threshold, 0, mask)
    
    # Apply mask
    enhanced = image + amount * mask
    
    return np.clip(enhanced, 0, 1)

def multi_scale_enhancement(image, scales=[1, 2, 4, 8]):
    """
    Multi-scale enhancement to bring out structures at different sizes
    """
    enhanced = image.copy()
    
    for scale in scales:
        # Create scale-specific enhancement
        blurred = ndimage.gaussian_filter(image, scale)
        detail = image - blurred
        enhanced += 0.3 * detail
    
    return np.clip(enhanced, 0, 1)

def noise_reduction(image, method='bilateral'):
    """
    Reduce noise while preserving structure
    """
    if method == 'bilateral':
        # Bilateral filter preserves edges
        return restoration.denoise_bilateral(image, sigma_color=0.1, sigma_spatial=1)
    elif method == 'wavelet':
        # Wavelet denoising
        return restoration.denoise_wavelet(image, method='BayesShrink', mode='soft')
    else:
        # Gaussian smoothing
        return ndimage.gaussian_filter(image, 0.5)

def create_publication_rgb(r_data, g_data, b_data, 
                          stretch_method='lupton',
                          enhance_detail=True,
                          reduce_noise=True,
                          color_balance=True):
    """
    Create publication-quality RGB composite
    """
    
    # Process each channel
    channels = []
    for data in [r_data, g_data, b_data]:
        if data is None:
            channels.append(np.zeros_like(r_data if r_data is not None else g_data))
            continue
            
        # Initial stretch
        stretched = professional_stretch(data, stretch_method, Q=10, alpha=1e-9)
        
        # Noise reduction
        if reduce_noise:
            stretched = noise_reduction(stretched, 'bilateral')
        
        # Detail enhancement
        if enhance_detail:
            stretched = unsharp_mask_pro(stretched, radius=1.5, amount=1.5)
            stretched = multi_scale_enhancement(stretched, scales=[1, 2, 4])
        
        channels.append(stretched)
    
    # Combine into RGB
    rgb = np.dstack(channels)
    
    # Color balance
    if color_balance:
        # Adjust each channel to have similar median values
        for i in range(3):
            channel = rgb[:,:,i]
            if np.max(channel) > 0:
                # Normalize to [0,1] range
                rgb[:,:,i] = channel / np.percentile(channel, 99.0)
    
    # Final saturation boost
    rgb = np.clip(rgb * 1.2, 0, 1)
    
    return rgb

def create_false_color_masterpiece(data, colormap='custom', enhance_structure=True):
    """
    Create false-color single-band image like Hubble/JWST releases
    """
    if data is None:
        return None
    
    # Professional stretch
    stretched = professional_stretch(data, 'lupton', Q=8, alpha=1e-10)
    
    if enhance_structure:
        # Multi-level enhancement
        stretched = unsharp_mask_pro(stretched, radius=2.0, amount=2.5)
        stretched = multi_scale_enhancement(stretched, scales=[1, 2, 4, 8])
        
        # Noise reduction
        stretched = noise_reduction(stretched, 'bilateral')
    
    # Apply colormap
    if colormap == 'custom':
        # Create custom dramatic colormap
        colors = ['#000000', '#1a0033', '#4d0066', '#cc0066', '#ff6600', '#ffcc00', '#ffffff']
        n_bins = 256
        cmap = LinearSegmentedColormap.from_list('custom', colors, N=n_bins)
    elif colormap == 'hubble_gold':
        colors = ['#000000', '#1a0f00', '#4d2b00', '#cc6600', '#ffaa00', '#ffdd99', '#ffffff']
        cmap = LinearSegmentedColormap.from_list('hubble_gold', colors, N=256)
    elif colormap == 'jwst_red':
        colors = ['#000000', '#1a0000', '#4d0000', '#990000', '#ff3300', '#ff9966', '#ffffff']
        cmap = LinearSegmentedColormap.from_list('jwst_red', colors, N=256)
    else:
        cmap = plt.cm.get_cmap(colormap)
    
    return cmap(stretched)

def create_news_worthy_image(galaxy_data, galaxy_name, output_style='dramatic'):
    """
    Create publication-quality images like those in astronomy news
    """
    
    filters = list(galaxy_data.keys())
    print(f"Available filters for {galaxy_name}: {filters}")
    
    # Load all available images
    images = {}
    for filt in filters:
        for file_type, file_path in galaxy_data[filt].items():
            if 'drc' in file_type or 'sci' in file_type:
                data, _ = load_fits_image(file_path)
                if data is not None:
                    images[filt] = data
                    print(f"Loaded {filt}: {data.shape}")
                break
    
    if not images:
        print("No images could be loaded!")
        return
    
    fig = plt.figure(figsize=(20, 15))
    
    if len(images) >= 3:
        # Multi-filter RGB composite
        filter_names = list(images.keys())[:3]
        r_img = images[filter_names[0]]
        g_img = images[filter_names[1]] 
        b_img = images[filter_names[2]]
        
        print(f"Creating RGB composite: R={filter_names[0]}, G={filter_names[1]}, B={filter_names[2]}")
        
        if output_style == 'dramatic':
            rgb = create_publication_rgb(r_img, g_img, b_img, 
                                       stretch_method='lupton',
                                       enhance_detail=True,
                                       reduce_noise=True)
        else:
            rgb = create_publication_rgb(r_img, g_img, b_img,
                                       stretch_method='adaptive',
                                       enhance_detail=False,
                                       reduce_noise=False)
        
        # Main large image
        ax1 = plt.subplot(2, 3, (1, 4))
        ax1.imshow(rgb, origin='lower')
        ax1.set_title(f'{galaxy_name.upper()}\nMulti-Wavelength Composite', 
                     fontsize=20, fontweight='bold', color='white')
        ax1.axis('off')
        
        # Individual filter panels
        for i, (filt, data) in enumerate(list(images.items())[:3]):
            ax = plt.subplot(2, 3, i + 4)
            false_color = create_false_color_masterpiece(data, 
                                                       colormap=['hubble_gold', 'jwst_red', 'custom'][i],
                                                       enhance_structure=True)
            ax.imshow(false_color, origin='lower')
            ax.set_title(f'{filt.upper()}', fontsize=14, fontweight='bold', color='white')
            ax.axis('off')
    
    else:
        # Single filter dramatic enhancement
        filt, data = next(iter(images.items()))
        
        # Create multiple dramatic versions
        styles = [
            ('Lupton Stretch', 'lupton', 'hubble_gold'),
            ('Adaptive Enhancement', 'adaptive', 'jwst_red'), 
            ('Multi-Scale Detail', 'power', 'custom')
        ]
        
        for i, (title, method, cmap) in enumerate(styles):
            ax = plt.subplot(1, 3, i + 1)
            
            if method == 'multi_scale':
                enhanced = professional_stretch(data, 'lupton')
                enhanced = multi_scale_enhancement(enhanced, scales=[1, 2, 4, 8])
                enhanced = unsharp_mask_pro(enhanced, radius=2.0, amount=3.0)
                false_color = create_false_color_masterpiece(enhanced, colormap=cmap, enhance_structure=False)
            else:
                false_color = create_false_color_masterpiece(data, colormap=cmap, enhance_structure=True)
            
            ax.imshow(false_color, origin='lower')
            ax.set_title(f'{galaxy_name.upper()}\n{title}', 
                        fontsize=16, fontweight='bold', color='white')
            ax.axis('off')
    
    # Set black background
    fig.patch.set_facecolor('black')
    plt.tight_layout()
    plt.subplots_adjust(wspace=0.02, hspace=0.1)
    
    # Save high-resolution version
    output_file = f"{galaxy_name}_publication_quality.png"
    plt.savefig(output_file, dpi=300, facecolor='black', bbox_inches='tight', 
                edgecolor='none', pad_inches=0.1)
    print(f"\n🎨 Saved publication-quality image: {output_file}")
    
    plt.show()
    
    return fig

def find_phangs_files(directory):
    """Find and categorize PHANGS FITS files"""
    pattern = os.path.join(directory, "hlsp_phangs*.fits")
    files = glob.glob(pattern)
    
    if not files:
        return {}
    
    galaxy_data = {}
    
    for file_path in files:
        filename = os.path.basename(file_path)
        parts = filename.split('_')
        
        try:
            galaxy = parts[4]  # e.g., ngc4535
            filter_name = parts[5]  # e.g., f275w
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
    """Create publication-quality galaxy images"""
    if len(sys.argv) != 2:
        print("Usage: python professional_imager.py /path/to/fits/directory")
        sys.exit(1)
    
    data_directory = sys.argv[1]
    
    print("🎨 Professional Galaxy Image Creator")
    print("="*50)
    print(f"Searching: {data_directory}")
    
    galaxy_data = find_phangs_files(data_directory)
    
    if not galaxy_data:
        print("No PHANGS data found!")
        sys.exit(1)
    
    for galaxy_name, galaxy_filters in galaxy_data.items():
        print(f"\n🌌 Creating publication images for {galaxy_name.upper()}")
        
        # Create dramatic version
        print("   → Creating dramatic enhancement...")
        create_news_worthy_image(galaxy_filters, galaxy_name, 'dramatic')
        
        # Create scientific version
        print("   → Creating scientific version...")  
        create_news_worthy_image(galaxy_filters, galaxy_name, 'scientific')
    
    print("\n✅ Publication-quality images complete!")
    print("Check the saved .png files for high-resolution versions!")

if __name__ == "__main__":
    main()