# PHANGS Galaxy Analysis Project

A Python toolkit for creating stunning, publication-quality images from PHANGS (Physics at High Angular resolution in Nearby GalaxieS) HST survey data. This project focuses on multi-wavelength analysis and visualization of barred spiral galaxies, with specialized tools for the NGC 4535 dataset.

## 🌌 Project Overview

This toolkit processes PHANGS FITS files to create professional astronomical images comparable to those seen in astronomy news articles and research publications. The project evolved from basic FITS viewing to sophisticated multi-wavelength composites that reveal different aspects of galaxy structure and star formation.

**Key Features:**

- 🔭 Multi-wavelength RGB composite creation
- 🎨 Professional image enhancement and stretching
- 📊 Automated galaxy structure analysis
- 💾 High-resolution publication-ready outputs
- 🌟 Specialized tools for UV-V-I and B-V-I combinations

## 🚀 Quick Start

### Best Results (Recommended)

For the highest quality images, use the UV-V-I focused tool:

```bash
python perfect_uv_v_i_creator.py /path/to/your/fits/files
```

This creates stunning composites highlighting star formation regions against overall galaxy structure.

## 📁 Project Structure

### Core Analysis Tools

| Script                      | Purpose                                  | Best For                         |
| --------------------------- | ---------------------------------------- | -------------------------------- |
| `phangs_analyzer.py`        | Basic FITS analysis and overview         | First-time data exploration      |
| `simple_galaxy_viewer.py`   | Clean, reliable visualization            | Getting familiar with your data  |
| `perfect_uv_v_i_creator.py` | ⭐ **Recommended** - UV-V-I masterpieces | Final publication-quality images |

### Specialized Tools

| Script                        | Purpose                         | Use Case                                |
| ----------------------------- | ------------------------------- | --------------------------------------- |
| `improved_galaxy_composer.py` | Multiple filter combinations    | Comparing different wavelength pairings |
| `focus_bvi_creator.py`        | B-V-I focused processing        | Classical optical combinations          |
| `perfect_galaxy_combos.py`    | Multiple combination comparison | Side-by-side analysis                   |

## 🛠️ Installation

### Requirements

```bash
pip install numpy matplotlib astropy scipy scikit-image
```

### Dependencies

- **numpy** - Array processing and mathematics
- **matplotlib** - Plotting and image creation
- **astropy** - FITS file handling and astronomical utilities
- **scipy** - Image processing and filtering
- **scikit-image** - Advanced image enhancement

## 📊 Expected Data Format

The tools expect PHANGS FITS files with naming convention:

```
hlsp_phangs-hst_hst_wfc3-uvis_ngc4535_f275w_v1_exp-drc-sci.fits
```

**Required File Types:**

- `exp-drc-sci` - Science data (primary)
- `drc` - Drizzle-corrected images (alternative)

**Supported Filters:**

- **F275W** - Near-UV (275 nm)
- **F336W** - U-band (336 nm)
- **F438W** - B-band (438 nm)
- **F555W** - V-band (555 nm)
- **F814W** - I-band (814 nm)

## 🎨 Output Examples

### UV-V-I Composite (Recommended)

```bash
python perfect_uv_v_i_creator.py ngc4535
```

**Creates:**

- `ngc4535_UV_V_I_clean_perfected.png` - Your preferred clean style
- `ngc4535_UV_V_I_clean_LARGE.png` - Extra high-resolution version
- `ngc4535_UV_V_I_comparison.png` - Multiple processing styles

**Scientific Value:**

- **Red (UV)**: Star formation regions and hot young stars
- **Green (V)**: Overall stellar distribution and galaxy structure
- **Blue (I)**: Older stellar populations and dust penetration

### Alternative Combinations

- **B-V-I**: Classical optical imaging for stellar populations
- **U-B-V**: Traditional color photography equivalent
- **UV-U-B**: Star formation emphasis

## 🔬 Scientific Background

### Why These Filter Combinations Work

**UV-V-I Combination:**

- Spans ~540 nm wavelength range (275-814 nm)
- Maximum contrast between star formation and stellar structure
- Excellent for studying spiral galaxy evolution
- Relates to recent discoveries of ancient barred spirals

**B-V-I Combination:**

- Traditional optical astronomy standard
- Optimal for stellar population analysis
- Good dust lane visibility
- Comparable to classical photographic plates

### Connection to Recent Research

This project was inspired by recent discoveries of ancient barred spiral galaxies:

- **J0107a** - 11.1 billion year old barred spiral discovered by ALMA
- **Event Horizon Telescope** - Multi-wavelength black hole imaging
- **JWST discoveries** - Early galaxy structure formation

## 📖 Usage Guide

### Basic Workflow

1. **Start Simple:**

   ```bash
   python simple_galaxy_viewer.py /path/to/fits/directory
   ```

2. **Explore Options:**

   ```bash
   python improved_galaxy_composer.py /path/to/fits/directory
   ```

3. **Create Masterpiece:**
   ```bash
   python perfect_uv_v_i_creator.py /path/to/fits/directory
   ```

### Command Line Arguments

All scripts take a single argument - the directory containing FITS files:

```bash
python script_name.py /path/to/fits/files
```

### Output Files

**High-Resolution Images:**

- 300 DPI PNG format
- Black backgrounds (professional astronomy style)
- Multiple size options (standard and LARGE versions)

**Comparison Charts:**

- Multiple processing methods side-by-side
- White backgrounds for presentations
- Lower resolution for quick viewing

## ⚙️ Advanced Configuration

### Customizing Processing

The tools include several adjustable parameters:

**Stretching Methods:**

- `zscale` - Percentile-based (recommended for clean results)
- `asinh` - Hyperbolic arcsine (handles wide dynamic range)
- `gamma` - Power law (good for enhancement)

**Color Balance:**

- `auto` - Statistical balancing
- `manual` - Preset weights for each filter

**Enhancement Options:**

- `unsharp_strength` - Detail enhancement (0.0-0.5 recommended)
- `saturation` - Color intensity (1.0-1.5 range)
- `brightness` - Overall brightness (0.8-1.2 range)

## 🐛 Troubleshooting

### Common Issues

**"No PHANGS FITS files found"**

- Check file naming convention
- Ensure files end in `.fits`
- Verify directory path

**"No science images found"**

- Look for files containing `drc` or `sci` in filename
- Check FITS file extensions (should be 1 or 0)

**Display errors/Tkinter issues**

- Tools save files regardless of display issues
- Check for saved PNG files in working directory
- Use `perfect_uv_v_i_creator.py` which avoids display problems

**Memory issues with large files**

- 7000x7000 pixel images require ~400MB RAM per filter
- Close other applications if needed
- Tools automatically manage memory by closing figures

### Performance Tips

- Start with `simple_galaxy_viewer.py` to verify data loads correctly
- Use `perfect_uv_v_i_creator.py` for final high-quality outputs
- Multiple filters will take 2-5 minutes to process
- High-resolution saves may take 30-60 seconds each

## 🌟 Best Practices

### For Publication-Quality Results

1. **Use clean processing:** Start with clean, conservative settings
2. **Focus on UV-V-I:** This combination provides the best scientific story
3. **Save multiple versions:** Keep both standard and LARGE outputs
4. **Check individual filters:** Verify each filter loaded correctly
5. **Document your process:** Note which settings produced best results

### Scientific Considerations

- **UV data is noisier:** Tools automatically handle this with specialized stretching
- **Color balance matters:** Automatic balancing usually works best
- **Structure vs. aesthetics:** Clean versions are better for scientific analysis
- **File sizes:** High-resolution outputs can be 50-100MB each

## 📚 References and Data Sources

**Data Source:**

- PHANGS-HST Treasury Survey
- Hubble Space Telescope observations
- Available through MAST (Mikulski Archive for Space Telescopes)

**Inspired by Recent Research:**

- J0107a ancient barred spiral discovery (ALMA 2025)
- Event Horizon Telescope multi-frequency imaging
- JWST early galaxy structure observations

**Processing Techniques:**

- Lupton et al. (2004) - Astronomical image stretching
- Professional astronomical image processing standards
- Publication-quality enhancement methods

## 🤝 Contributing

This project evolved through experimentation with different processing approaches. Feel free to:

- Modify stretching parameters for your specific data
- Add new filter combinations
- Improve error handling for different FITS formats
- Add support for other survey data

## 📄 License

Open source - feel free to use and modify for your astronomical research and visualization projects.

---

**Created for exploring the beautiful structure of NGC 4535 and other PHANGS survey galaxies. Happy imaging! 🌌**
