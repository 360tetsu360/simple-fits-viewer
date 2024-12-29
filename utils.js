function parseFITSImage(arrayBuffer, dataView) {
    
    console.time("parseFITSImage");

    // Very basic FITS header parsing
    let headerText = "";
    let offset = 0;
    const headerSize = 2880;
    while (true) {
        const block = new TextDecoder().decode(
            arrayBuffer.slice(offset, offset + headerSize)
        );
        headerText += block;
        offset += headerSize;
        if (block.trim().endsWith("END")) break;
    }

    // Parse Header Keywords
    const headerLines = headerText.match(/.{1,80}/g); // Split into 80-char lines
    const header = {};
    for (const line of headerLines) {
        const keyword = line.substring(0, 8).trim();
        const value = line.substring(10, 80).trim();
        if (keyword === "END") break;
        header[keyword] = value;
    }
    console.timeLog("parseFITSImage", "parseFITSHeader");

    const width = parseInt(header["NAXIS1"], 10);
    const height = parseInt(header["NAXIS2"], 10);
    const bitpix = parseInt(header["BITPIX"], 10);
    const bscale = parseFloat(header["BSCALE"]) || 1;
    const bzero = parseFloat(header["BZERO"]) || 0;

    // Parse Image Data
    const dataSize = width * height;
    const bytesPerPixel = Math.abs(bitpix) / 8;

    // Use a typed array for image data
    let data;
    if (bitpix === 8 || bitpix === 16 || bitpix === 32) {
        data = new Int32Array(dataSize);
    } else if (bitpix === -32) {
        data = new Float32Array(dataSize);
    } else if (bitpix === -64) {
        data = new Float64Array(dataSize);
    } else {
        throw new Error(`Unsupported BITPIX: ${bitpix}`);
    }

    for (let i = 0; i < dataSize; i++) {
        if (bitpix === 8) {
            data[i] = dataView.getUint8(offset) * bscale + bzero;
        } else if (bitpix === 16) {
            data[i] = dataView.getInt16(offset, false) * bscale + bzero;
        } else if (bitpix === 32) {
            data[i] = dataView.getInt32(offset, false) * bscale + bzero;
        } else if (bitpix === -32) {
            data[i] = dataView.getFloat32(offset, false) * bscale + bzero;
        } else if (bitpix === -64) {
            data[i] = dataView.getFloat64(offset, false) * bscale + bzero;
        }
        offset += bytesPerPixel;
    }
    console.timeLog("parseFITSImage", "parseFITSImageData");

    // Normalize Data for Display
    const { vmin, vmax } = zscale(data);
    console.timeLog("parseFITSImage", "zscale");
    // const normalizedData = data.map(
    //     (value) => ((value - vmin) / (vmax - vmin)) * 255
    // );
    const scale = 255 / (vmax - vmin);
    const _offset = -vmin * scale;
    const normalizedData = new Array(data.length);

    for (let i = 0; i < data.length; i++) {
        normalizedData[i] = data[i] * scale + _offset;
    }
    console.timeLog("parseFITSImage", "normalizeData");

    console.timeEnd("parseFITSImage", "parseFITSImage done");

    // console.log(header, normalizedData);
    return [header, normalizedData, width, height, data];
}


function zscale(
    values,
    n_samples = 1000,
    contrast = 0.25,
    max_reject = 0.5,
    min_npixels = 5,
    krej = 2.5,
    max_iterations = 5
) {
    console.time("zscale");

    // Sample the image
    const stride = Math.max(1, Math.floor(values.length / n_samples));
    const samples = [];
    for (let i = 0; i < values.length && samples.length < n_samples; i += stride) {
        samples.push(values[i]);
    }
    console.timeLog("zscale", "sampleImage");

    // Sort in-place to avoid extra memory usage
    samples.sort((a, b) => a - b);
    console.timeLog("zscale", "sortSamples");

    const npix = samples.length;
    let vmin = samples[0];
    let vmax = samples[npix - 1];

    // Precompute x values
    const x = new Array(npix);
    for (let i = 0; i < npix; i++) {
        x[i] = i;
    }
    console.timeLog("zscale", "precomputeX");

    let ngoodpix = npix;
    let last_ngoodpix = ngoodpix + 1;

    // Initialize bad pixels mask
    const badpix = new Array(npix).fill(false);

    const minpix = Math.max(min_npixels, Math.floor(npix * max_reject));
    let fit = { slope: 0, intercept: 0 };
    console.timeLog("zscale", "initializeBadPixelsMask");

    for (let iter = 0; iter < max_iterations; iter++) {
        if (ngoodpix >= last_ngoodpix || ngoodpix < minpix) break;

        fit = linearFit(x, samples, badpix);
        // Compute fitted values and residuals using loops
        const fitted = new Array(npix);
        const flat = new Array(npix);
        for (let i = 0; i < npix; i++) {
            fitted[i] = fit.slope * x[i] + fit.intercept;
            flat[i] = samples[i] - fitted[i];
        }

        // Compute threshold for k-sigma clipping
        const goodPixels = [];
        for (let i = 0; i < npix; i++) {
            if (!badpix[i]) goodPixels.push(flat[i]);
        }
        const sigma = std(goodPixels);
        const threshold = krej * sigma;

        // Update badpix mask
        ngoodpix = 0;
        for (let i = 0; i < npix; i++) {
            if (Math.abs(flat[i]) > threshold) {
                badpix[i] = true;
            } else {
                badpix[i] = false;
                ngoodpix++;
            }
        }

        last_ngoodpix = ngoodpix;
    }
    console.timeLog("zscale", "kSigmaClipping");

    if (ngoodpix >= minpix) {
        let slope = fit.slope;
        if (contrast > 0) {
            slope = slope / contrast;
        }
        const center_pixel = Math.floor((npix - 1) / 2);
        const median = medianValue(samples);
        vmin = Math.max(vmin, median - (center_pixel - 1) * slope);
        vmax = Math.min(vmax, median + (npix - center_pixel) * slope);
    }
    console.timeLog("zscale", "updateMinMax");

    return { vmin, vmax };
}

function linearFit(x, y, badpix) {
    // Optimized linear fit using loops
    let sumX = 0,
        sumY = 0,
        sumXY = 0,
        sumX2 = 0,
        n = 0;
    for (let i = 0; i < x.length; i++) {
        if (!badpix[i]) {
            const xi = x[i];
            const yi = y[i];
            sumX += xi;
            sumY += yi;
            sumXY += xi * yi;
            sumX2 += xi * xi;
            n++;
        }
    }
    const denominator = n * sumX2 - sumX * sumX;
    const slope = (n * sumXY - sumX * sumY) / denominator;
    const intercept = (sumY - slope * sumX) / n;
    return { slope, intercept };
}

function std(arr) {
    // Optimized standard deviation calculation
    let mean = 0;
    for (let i = 0; i < arr.length; i++) {
        mean += arr[i];
    }
    mean /= arr.length;
    let variance = 0;
    for (let i = 0; i < arr.length; i++) {
        const diff = arr[i] - mean;
        variance += diff * diff;
    }
    variance /= arr.length;
    return Math.sqrt(variance);
}

function medianValue(arr) {
    // Optimized median calculation using Quickselect algorithm
    const n = arr.length;
    const k = Math.floor(n / 2);
    return quickSelect(arr, k);
}

function quickSelect(arr, k) {
    // In-place Quickselect algorithm
    let left = 0;
    let right = arr.length - 1;
    while (left <= right) {
        const pivotIndex = partition(arr, left, right);
        if (pivotIndex === k) {
            return arr[k];
        } else if (pivotIndex < k) {
            left = pivotIndex + 1;
        } else {
            right = pivotIndex - 1;
        }
    }
}

function partition(arr, left, right) {
    const pivotValue = arr[right];
    let pivotIndex = left;
    for (let i = left; i < right; i++) {
        if (arr[i] < pivotValue) {
            [arr[i], arr[pivotIndex]] = [arr[pivotIndex], arr[i]];
            pivotIndex++;
        }
    }
    [arr[right], arr[pivotIndex]] = [arr[pivotIndex], arr[right]];
    return pivotIndex;
}

function convolve(arr, kernel) {
    // Optimized convolution using loops
    const result = new Array(arr.length).fill(false);
    const kernelLength = kernel.length;
    for (let i = 0; i < arr.length; i++) {
        if (arr[i]) {
            for (let j = 0; j < kernelLength; j++) {
                const idx = i + j;
                if (idx < arr.length) {
                    result[idx] = true;
                }
            }
        }
    }
    return result;
}

function formatNumber(num, precision) {
    if (Math.floor(num) === num) {
        return num; // return as is, when it's an integer
    } else {
        return num.toFixed(precision); // use toFixed when there are decimals
    }
}

const deg2arcsec = (deg) => deg * 3600.0;
const deg2rad = (deg) => deg * (Math.PI / 180);
const radec2x = (r,d) => Math.cos(d)*Math.cos(r);
const radec2y = (r,d) => Math.cos(d)*Math.sin(r);
const radec2z = (r,d) => Math.sin(d);
const radec2xyz = (r,d) => [radec2x(r,d), radec2y(r,d), radec2z(r,d)];
const radecdeg2xyz = (r,d) => radec2xyz(deg2rad(r), deg2rad(d));

function normalize(x, y, z) {
	let invl = 1.0 / Math.sqrt(x*x + y*y + z*z);
    return [x * invl, y * invl, z * invl];
}

const rad2deg = (rad) => rad * (180 / Math.PI);
const z2dec = (z) => Math.asin(z);
function xy2ra(x, y) {
    let a = Math.atan2(y, x);
    if (a < 0) {
        a += 2.0 * Math.PI;
    }
    return a;
}

function xyz2radec(x, y, z) {
    return [xy2ra(x, y), z2dec(z)];
}

function xyz2radecdeg(xyz) {
    let [ra, dec] = xyz2radec(xyz[0], xyz[1], xyz[2]);
    return [rad2deg(ra), rad2deg(dec)];
}

function wcs_pixel_center_for_size(size) {
    return 0.5 + 0.5 * size;
}

function star_coords(s, r, tangent) {
    let x = 0;
    let y = 0;
    // As used by the sip.c code, this does the TAN projection
    // (if "tangent" is TRUE; SIN projection otherwise)
    // r: CRVAL
    // s: RA,Dec to be projected
    // ASSUME r,s are unit vectors
    // sdotr:  s dot r = |r||s| cos(theta) = cos(theta)
    let sdotr = s[0] * r[0] + s[1] * r[1] + s[2] * r[2];
    if (sdotr <= 0.0) {
        // on the opposite side of the sky
        return null;
    }
    if (r[2] == 1.0) {
        // North pole
        let inv_s2 = 1.0 / s[2];
        if (tangent) {
            x = s[0] * inv_s2;
            y = s[1] * inv_s2;
        } else {
            x = s[0];
            y = s[1];
        }
    } else if (r[2] == -1.0) {
        // South pole
        let inv_s2 = 1.0 / s[2];
        if (tangent) {
            x = -s[0] * inv_s2;
            y =  s[1] * inv_s2;
        } else {
            x = -s[0];
            y =  s[1];
        }
    } else {
        // eta is a vector perpendicular to r pointing in the direction
        // of increasing RA.  eta_z = 0 by definition.
        let etax = -r[1];
        let etay =  r[0];
        let eta_norm = Math.hypot(etax, etay);
        let inv_en = 1.0 / eta_norm;
        etax *= inv_en;
        etay *= inv_en;

        // xi =  r cross eta, a vector pointing northwards,
        // in direction of increasing DEC
        let xix = -r[2] * etay;
        let xiy =  r[2] * etax;
        let xiz =  r[0] * etay - r[1] * etax;

        // project s-r onto eta and xi.  No need to subtract r from s, though,
        // since eta and xi are orthogonal to r by construction.
        x = (s[0] * etax + s[1] * etay             );
        y = (s[0] *  xix + s[1] *  xiy + s[2] * xiz);

        // The "inv_sdotr" applies the TAN scaling
        if (tangent) {
            let inv_sdotr = 1.0 / sdotr;
            x *= inv_sdotr;
            y *= inv_sdotr;
        }
    }
    return [x, y];
}

function invert_2by2_arr(a) {
    let ainv = new Array(2).fill(0).map(() => new Array(2).fill(0));
    let det = a[0][0] * a[1][1] - a[0][1] * a[1][0];

    if (det == 0.0) {
        return null;
    }

    let inv_det = 1.0 / det;
    ainv[0][0] =  a[1][1] * inv_det;
    ainv[0][1] = -a[0][1] * inv_det;
    ainv[1][0] = -a[1][0] * inv_det;
    ainv[1][1] =  a[0][0] * inv_det;
    return ainv;
}


function parseWCSPolynomial(header, name, order) {
    const data = Array.from({ length: 10 }, () => Array(10).fill(0));

    for (let i = 0; i < order; i++) {
        for (let j = 0; j < order; j++) {
            const key = `${name}_${i}_${j}`;
            if (header[key] !== undefined) {
                data[i][j] = parseFloat(header[key].split("/")[0]);
            }
        }
    }

    return data;
}

class WCS {
    // Parse WCS header
    // check if wcs is included
    constructor(header) {
        this.wcsaxes = parseInt(header["WCSAXES"].split("/")[0], 10);
        this.ctype1 = header["CTYPE1"].split("/")[0];
        this.ctype2 = header["CTYPE2"].split("/")[0];

        if (this.ctype1 == "RA---SIN-SIP" && this.ctype2 == "DEC---SIN-SIP") {
            this.sin = true;
        } else if (this.ctype1 == "RA---TAN-SIP" && this.ctype2 == "DEC---TAN-SIP") {
            this.sin = false;
        } else {
            console.error(`Unsupported wcs format: ${this.ctype1}`);
        }

        this.equinox = parseFloat(header["EQUINOX"].split("/")[0]);
        this.lonpole = parseFloat(header["LONPOLE"].split("/")[0]);
        this.latpole = parseFloat(header["LATPOLE"].split("/")[0]);

        this.crval = [];
        this.crval.push(parseFloat(header["CRVAL1"].split("/")[0]));
        this.crval.push(parseFloat(header["CRVAL2"].split("/")[0]));

        this.crpix = [];
        this.crpix.push(parseFloat(header["CRPIX1"].split("/")[0]));
        this.crpix.push(parseFloat(header["CRPIX2"].split("/")[0]));

        this.cunit1 = header["CUNIT1"].split("/")[0];
        this.cunit2 = header["CUNIT2"].split("/")[0];

        this.cd = new Array(2).fill(0).map(() => new Array(2).fill(0));
        this.cd[0][0] = parseFloat(header["CD1_1"].split("/")[0]);
        this.cd[0][1] = parseFloat(header["CD1_2"].split("/")[0]);
        this.cd[1][0] = parseFloat(header["CD2_1"].split("/")[0]);
        this.cd[1][1] = parseFloat(header["CD2_2"].split("/")[0]);

        this.imagew = parseInt(header["IMAGEW"].split("/")[0], 10);
        this.imageh = parseInt(header["IMAGEH"].split("/")[0], 10);

        this.a_order = parseInt(header["A_ORDER"].split("/")[0], 10);
        this.a = parseWCSPolynomial(header, "A", this.a_order);
        this.b_order = parseInt(header["B_ORDER"].split("/")[0], 10);
        this.b = parseWCSPolynomial(header, "B", this.b_order);
        this.ap_order = parseInt(header["AP_ORDER"].split("/")[0], 10);
        this.ap = parseWCSPolynomial(header, "AP", this.ap_order);
        this.bp_order = parseInt(header["BP_ORDER"].split("/")[0], 10);
        this.bp = parseWCSPolynomial(header, "BP", this.ap_order);
    }

    sip_get_radec_bounds(stepsize) {
        let [rac, decc] = this.sip_get_radec_center();

        let [ramin, ramax, decmin, decmax] = this.sip_walk_image_boundary(stepsize, rac, decc);

        // Check for poles...
        // north pole
        if (this.sip_is_inside_image(0, 90)) {
            ramin = 0;
            ramax = 360;
            decmax = 90;
        }
        if (this.sip_is_inside_image(0, -90)) {
            ramin = 0;
            ramax = 360;
            decmin = -90;
        }

        return [ramin, ramax, decmin, decmax];
    }

    sip_is_inside_image(ra, dec) {
        let xy = this.sip_radec2pixelxy(ra, dec);
        if (xy == null) {
            return false;
        }
        return this.tan_pixel_is_inside_image(xy[0], xy[1]);
    }

    sip_radec2pixelxy(ra, dec) {
        let xy = this.tan_radec2pixelxy(ra, dec);
        if (xy == null) {
            return null;
        }
        return this.sip_pixel_undistortion(xy[0], xy[1]);
    }

    sip_walk_image_boundary(stepsize, rac, decc) {
        // Walk the perimeter of the image in steps of stepsize pixels
        let ramin = rac;
        let ramax = rac;
        let decmin = decc;
        let decmax = decc;

        let w = this.imagew;
        let h = this.imageh;
        
        let xmin = 0.5;
        let xmax = w + 0.5;
        let ymin = 0.5;
        let ymax = h + 0.5;
        let offsetx = [xmin, xmax, xmax, xmin];
        let offsety = [ymin, ymin, ymax, ymax];
        let stepx = [+stepsize, 0, -stepsize, 0];
        let stepy = [0, +stepsize, 0, -stepsize];
        let nsteps = [Math.ceil(w/stepsize), Math.ceil(h/stepsize), Math.ceil(w/stepsize), Math.ceil(h/stepsize) ];

        for (let side = 0; side < 4; side++) {
            for (let i = 0; i < nsteps[side]; i++) {
                let x = Math.min(xmax, Math.max(xmin, offsetx[side] + i * stepx[side]));
                let y = Math.min(ymax, Math.max(ymin, offsety[side] + i * stepy[side]));

                let [ra, dec] = this.sip_pixelxy2radec(x, y);
                decmin = Math.min(decmin, dec);
                decmax = Math.max(decmax, dec);
                if (ra - rac > 180) {
                    // wrap-around: racenter < 180, ra has gone < 0 but been wrapped around to > 180
                    ra -= 360;
                }
                if (rac - ra > 180) {
                    // wrap-around: racenter > 180, ra has gone > 360 but wrapped around to > 0.
                    ra += 360;
                }
            
                ramin = Math.min(ramin, ra);
                ramax = Math.max(ramax, ra);
            }
        }

        return [ramin, ramax, decmin, decmax];
    }

    sip_get_radec_center() {
        let px = wcs_pixel_center_for_size(this.imagew);
        let py = wcs_pixel_center_for_size(this.imageh);
        return this.sip_pixelxy2radec(px, py);
    }

    has_distortions() {
        return this.a_order >= 0;
    }

    sip_pixelxy2radec(px, py) {
        if (this.has_distortions()) {
            let [u, v] = this.sip_distortion(px, py);
            // Run a normal TAN conversion on the distorted pixel coords.
            return this.tan_pixelxy2radec(u, v);
        } else {
            // Run a normal TAN conversion
            return this.tan_pixelxy2radec(px, py);
        }
    }

    sip_distortion(px, py) {
        // Get pixel coordinates relative to reference pixel
        let u = px - this.crpix[0];
        let v = py - this.crpix[1];
        let xy = this.sip_calc_distortion(u, v);
        xy[0] += this.crpix[0];
        xy[1] += this.crpix[1];
        return xy;
    }

    sip_calc_distortion(u, v) {    
        let fuv = 0.0;
        let guv = 0.0;

        // avoid using pow() function
        const powu = new Array(10).fill(0.0);
        const powv = new Array(10).fill(0.0);

        powu[0] = 1.0;
        powu[1] = u; 
        powv[0] = 1.0;
        powv[1] = v; 

        for (let i = 2; i <= Math.max(this.a_order, this.b_order); i++) {
            powu[i] = powu[i - 1] * u; // u^i = u^(i-1) * u
            powv[i] = powv[i - 1] * v; // v^i = v^(i-1) * v
        }
        
        for (let i = 0; i <= this.a_order; i++) {
            for (let j = 0; j <= this.a_order; j++) {
                // We include all terms, even the constant and linear ones; the standard
                // isn't clear on whether these are allowed or not.
                if (i+j <= this.a_order) {
                    fuv += this.a[i][j] * powu[i] * powv[j];
                }
            }
        }

        for (let i = 0; i<=this.b_order; i++) {
            for (let j=0; j<=this.b_order; j++) {
                if (i+j <= this.b_order) {
                    guv += this.b[i][j] * powu[i] * powv[j];
                }
            }
        }

        return [u + fuv, v + guv];
    }

    tan_pixelxy2iwc(px, py) {
        // Get pixel coordinates relative to reference pixel
        let u = px - this.crpix[0];
        let v = py - this.crpix[1];

        // Get intermediate world coordinates
        let x = this.cd[0][0] * u + this.cd[0][1] * v;
        let y = this.cd[1][0] * u + this.cd[1][1] * v;

        return [x, y]
    }

    tan_iwc2xyzarr(x, y)
    {
        let ix,iy,norm;
        let jx,jy,jz;
        let xyz = [0, 0, 0];
    
        // Mysterious factor of -1 correcting for vector directions below.
        x = -deg2rad(x);
        y =  deg2rad(y);

        // Take r to be the threespace vector of crval
        let [rx, ry, rz] = radecdeg2xyz(this.crval[0], this.crval[1]);
        //printf("rx=%lf ry=%lf rz=%lf\n",rx,ry,rz);
    
        // FIXME -- what about *near* the poles?
        if (rx == 1.0) {
            // North pole
            ix = -1.0;
            iy = 0.0;
        } else if (rz == -1.0) {
            // South pole
            ix = -1.0;
            iy = 0.0;
        } else {
            // Form i = r cross north pole (0,0,1)
            ix = ry;
            iy = -rx;
            // iz = 0
            norm = Math.hypot(ix, iy);
            ix /= norm;
            iy /= norm;
            //printf("ix=%lf iy=%lf iz=0.0\n",ix,iy);
            //	printf("r.i = %lf\n",ix*rx+iy*ry);
        }
    
        // Form j = i cross r;   iz=0 so some terms drop out
        jx = iy * rz;
        jy =         - ix * rz;
        jz = ix * ry - iy * rx;
        // norm should already be 1, but normalize anyway
        let [jx_, jy_, jz_] = normalize(jx, jy, jz);
        //	printf("jx=%lf jy=%lf jz=%lf\n",jx,jy,jz);
        //	printf("r.j = %lf\n",jx*rx+jy*ry+jz*rz);
        //	printf("i.j = %lf\n",ix*jx+iy*jy);
    
        if (this.sin) {
            console.assert((x*x + y*y) < 1.0);
            // Figure out what factor of r we have to add in to make the resulting length = 1
            let rfrac = Math.sqrt(1.0 - (x*x + y*y));
            // Don't scale the projected x,y positions, just add in the right amount of r to
            // bring it onto the unit sphere
            xyz[0] = ix*x + jx_*y + rx * rfrac;
            xyz[1] = iy*x + jy_*y + ry * rfrac;
            xyz[2] =        jz_*y + rz * rfrac; // iz = 0
            return xyz;
        } else {
            // Form the point on the tangent plane relative to observation point,
            xyz[0] = ix*x + jx_*y + rx;
            xyz[1] = iy*x + jy_*y + ry;
            xyz[2] =        jz_*y + rz; // iz = 0
            // and normalize back onto the unit sphere
            return normalize(xyz[0], xyz[1], xyz[2]);
        }
    }
    //
    tan_pixelxy2xyzarr(px, py) {
        let [x, y] = this.tan_pixelxy2iwc(px, py);
        return this.tan_iwc2xyzarr(x, y);
    }

    tan_pixelxy2radec(px, py) {
        let xyz = this.tan_pixelxy2xyzarr(px, py);
        return xyz2radecdeg(xyz);
    }

    tan_radec2pixelxy(a, d) {
        let xyzpt = radecdeg2xyz(a,d);
        return this.tan_xyzarr2pixelxy(xyzpt);
    }

    tan_xyzarr2pixelxy(xyz) {
        let iw = this.tan_xyzarr2iwc(xyz);
        if (iw == null) {
            return null;
        }
        return this.tan_iwc2pixelxy(iw[0], iw[1]);
    }

    tan_xyzarr2iwc(xyz) {
        // FIXME be robust near the poles
        // Calculate intermediate world coordinates (x,y) on the tangent plane
        let xyzcrval = radecdeg2xyz(this.crval[0], this.crval[1]);

        let iw = star_coords(xyz, xyzcrval, !this.sin)
        if (iw == null) {
            return null;
        }

        let iwcx = rad2deg(iw[0]);
        let iwcy = rad2deg(iw[1]);
        return [iwcx, iwcy];
    }

    tan_iwc2pixelxy(x, y) {
        // Invert CD
        let cdi = invert_2by2_arr(this.cd);

        // Linear pixel coordinates
        let u = cdi[0][0]*x + cdi[0][1]*y;
        let v = cdi[1][0]*x + cdi[1][1]*y;

        // Re-add crpix to get pixel coordinates
        let px = u + this.crpix[0];
        let py = v + this.crpix[1];

        return [px, py];
    }

    tan_pixel_is_inside_image(x, y) {
        return (x >= 1 && x <= this.imagew && y >= 1 && y <= this.imageh);
    }
    

    sip_pixel_undistortion(x, y) {
        if (!this.has_distortions()) {
            return [x, y];
        }
        // Sanity check:
        if (this.a_order != 0 && this.ap_order == 0) {
            console.error("suspicious inversion; no inverse SIP coeffs yet there are forward SIP coeffs");
        }
    
        // Get pixel coordinates relative to reference pixel
        let u = x - this.crpix[0];
        let v = y - this.crpix[1];
        [x, y] = this.sip_calc_inv_distortion(u, v);
        x += this.crpix[0];
        y += this.crpix[1];
        return [x, y]; 
    }

    sip_calc_inv_distortion(u, v) {
        let fuv = 0.0;
        let guv = 0.0;

        // avoid using pow() function
        const powu = new Array(10).fill(0.0);
        const powv = new Array(10).fill(0.0);

        powu[0] = 1.0;
        powu[1] = u; 
        powv[0] = 1.0;
        powv[1] = v; 

        for (let i = 2; i <= Math.max(this.a_order, this.b_order); i++) {
            powu[i] = powu[i - 1] * u; // u^i = u^(i-1) * u
            powv[i] = powv[i - 1] * v; // v^i = v^(i-1) * v
        }

        for (let i = 0; i <= this.ap_order; i++) {
            for (let j = 0; j <= this.ap_order; j++) {
                if (i+j <= this.ap_order) {
                    fuv += this.ap[i][j] * powu[i] * powv[j];
                }
            }
        }

        return [u + fuv, v + guv];
    }

    tan_pixel_scale() {
        let scale = deg2arcsec(Math.sqrt(Math.abs(this.tan_det_cd())));
        return scale;
    }

    tan_det_cd() {
        return (this.cd[0][0]*this.cd[1][1] - this.cd[0][1]*this.cd[1][0]);
    }
}
