# Circle Size Legend Scaling Guide

## 📋 Overview
This guide helps you ensure perfect size fidelity when manually adding the circle size legend to your mega visual plots in external applications (Adobe Illustrator, Inkscape, etc.).

## 🎯 Files Generated
- `circle_size_legend.png/pdf` - The standalone legend to add to your figures
- `circle_size_reference.png/pdf` - Reference plot with known ratios for size comparison

## 📏 Scaling Fidelity Method

### Method 1: Reference Plot Comparison (Recommended)
1. **Open both files** in your graphics application:
   - Your main mega visual plot (16S or 18S)
   - The `circle_size_reference.png` file

2. **Compare circle sizes** for the labeled ratios:
   - Look for circles with 1× ratio (smallest)
   - Look for circles with 10× ratio (medium)  
   - Look for circles with 100× ratio (largest)

3. **Adjust legend scale** until the reference circles match your main plot circles exactly

4. **Apply same scaling** to the `circle_size_legend.png`

### Method 2: DPI-Based Scaling
1. **Check your main plot DPI**: Both scripts generate at 300 DPI
2. **Import legend at same DPI**: Use 300 DPI when importing the legend
3. **Import at 100% scale**: This should give perfect size matching
4. **Fine-tune if needed**: Minor adjustments may be needed due to application differences

### Method 3: Manual Measurement
1. **Measure a known circle** in your main plot (use ruler tool)
2. **Find same ratio** in the reference plot
3. **Measure that circle** and calculate scale factor
4. **Apply scale factor** to the legend

## 🎨 Legend Specifications

### Technical Details
- **Size Range**: 10-22 (same as main plots)
- **Transformation**: Square root of Genome/Isolate ratio
- **Scaling Formula**: `3 + 9 * (sqrt(ratio) - min) / (max - min)`
- **DPI**: 300 (matches main plots)
- **Background**: White
- **Dimensions**: 6" × 8"

### Representative Ratios Shown
- **1.0×**: Minimum meaningful ratio (smallest circle)
- **10.0×**: Medium ratio (medium circle)
- **100.0×**: High ratio (largest circle)

## 🔧 Application-Specific Tips

### Adobe Illustrator
1. Place both main plot and reference plot
2. Use the "Scale" tool with "Uniform" scaling
3. Hold Shift while scaling to maintain proportions
4. Use "Object > Transform > Scale" for precise percentage scaling

### Inkscape
1. Import at same DPI (300)
2. Use "Object > Transform" for precise scaling
3. Lock width/height ratio during scaling
4. Use "View > Zoom > 1:1" to check actual size

### Adobe Photoshop
1. Check "Image > Image Size" for both files
2. Ensure same resolution (300 pixels/inch)
3. Use "Edit > Transform > Scale" with Shift key
4. Check "Maintain Aspect Ratio" in options

## ✅ Verification Checklist

- [ ] Reference plot circles match main plot circles for same ratios
- [ ] Legend circles follow the same size progression
- [ ] Smallest circle (1×) is clearly visible but not too large
- [ ] Largest circle (100×) is prominent but not overwhelming
- [ ] Size progression looks smooth and logical
- [ ] Text labels are readable and properly positioned

## 🎯 Expected Results

When properly scaled, you should see:
- **Smooth size progression** from 1× to 100×
- **Visually intuitive scaling** (10× should look ~3× larger than 1×)
- **Perfect match** with circles in your main plots
- **Professional appearance** suitable for publication

## 🚨 Troubleshooting

### Circles Too Small
- Increase legend scale by 10-20%
- Check if main plot was scaled down during export

### Circles Too Large  
- Decrease legend scale by 10-20%
- Verify DPI settings match between files

### Sizes Don't Match Reference
- Re-export main plot at 300 DPI
- Ensure no compression was applied during export
- Check that both files use same units (inches/pixels)

### Text Overlaps Circles
- The legend is designed with adequate spacing
- If overlap occurs, the legend may be over-scaled
- Reduce scale slightly or increase legend canvas size

## 📊 Quality Control

The circle size legend uses the **exact same mathematical transformation** as both the 16S and 18S mega visual scripts:

```r
# Same transformation used in all scripts
Genome_Isolate_Ratio <- pmax(NCBI_Genome_Count / pmax(Isolate_Count, 1), 1)
Circle_Size_Raw <- sqrt(Genome_Isolate_Ratio)  
Circle_Size <- 3 + 9 * (Circle_Size_Raw - min) / (max - min)
scale_size_continuous(range = c(10, 22))
```

This ensures **perfect mathematical consistency** across all your visualizations.
