#
# LSST Data Management System
# Copyright 2008, 2009, 2010 LSST Corporation.
#
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <http://www.lsstcorp.org/LegalNotices/>.
#

import lsst.afw.image as afwImage
import lsst.meas.algorithms as measAlg
import lsst.afw.cameraGeom as cameraGeom
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from . import isr
from .isrLib import UnmaskedNanCounterF
from .assembleCcdTask import AssembleCcdTask

class IsrTaskConfig(pexConfig.Config):
    doWrite = pexConfig.Field(
        dtype = bool,
        doc = "Persist visitCCD?",
        default = True,
    )
    assembleCcd = pexConfig.ConfigurableField(
        target = AssembleCcdTask,
        doc = "Assemble CCD",
    )
    fwhm = pexConfig.Field(
        dtype = float,
        doc = "FWHM of PSF (arcsec)",
        default = 1.0,
    )
    #This is needed for both the detection and correction aspects
    saturatedMaskName = pexConfig.Field(
        dtype = str,
        doc = "Name of mask plane to use in saturation detection",
        default = "SAT",
    )
    flatScalingType = pexConfig.ChoiceField(
        dtype = str,
        doc = "The method for scaling the flat on the fly.",
        default = 'USER',
        allowed = {
            "USER":   "Scale by flatUserScale",
            "MEAN":   "Scale by the inverse of the mean",
            "MEDIAN": "Scale by the inverse of the median",
        },
    )
    flatUserScale = pexConfig.Field(
        dtype = float,
        doc = "If flatScalingType is 'USER' then scale flat by this amount; ignored otherwise",
        default = 1.0,
    )
    overscanFitType = pexConfig.ChoiceField(
        dtype = str,
        doc = "The method for fitting the overscan bias level.",
        default = 'MEDIAN',
        allowed = {
            "POLY": "Fit polynomial to the longest axis of the overscan region",
            "MEAN": "Correct using the mean of the overscan region",
            "MEDIAN": "Correct using the median of the overscan region",
        },
    )
    overscanPolyOrder = pexConfig.Field(
        dtype = int,
        doc = "Order of polynomial to fit if overscan fit type is POLY",
        default = 1,
    )
    growSaturationFootprintSize = pexConfig.Field(
        dtype = int,
        doc = "Number of pixels by which to grow the saturation footprints",
        default = 1,
    )
    growDefectFootprintSize = pexConfig.Field(
        dtype = int,
        doc = "Number of pixels by which to grow the defect (bad and nan) footprints",
        default = 1,
    )
    setGainAssembledCcd = pexConfig.Field(
        dtype = bool,
        doc = "update exposure metadata in the assembled ccd to reflect the effective gain of the assembled chip",
        default = True,
    )
    keysToRemoveFromAssembledCcd = pexConfig.ListField(
        dtype = str,
        doc = "fields to remove from the metadata of the assembled ccd.",
        default = [],
    )
    reNormAssembledCcd = pexConfig.Field(
        dtype = bool,
        doc = "renormalize the assembled chips to have unity gain.  False if setGain is False",
        default = True,
    )
    
    
class IsrTask(pipeBase.Task):
    ConfigClass = IsrTaskConfig
    _DefaultName = "isr"

    def __init__(self, *args, **kwargs):
        pipeBase.Task.__init__(self, *args, **kwargs)
        self.makeSubtask("assembleCcd")
        self.transposeForInterpolation = False

    @pipeBase.timeMethod
    def run(self, sensorRef):
        """Do instrument signature removal on an exposure
        
        Correct for saturation, bias, overscan, dark, flat..., perform CCD assembly,
        and interpolate over defects and saturated pixels.
        
        Persist the ISR-corrected exposure as "visitCcd" if config.doWrite is True.

        @param sensorRef daf.persistence.butlerSubset.ButlerDataRef of the data to be processed
        @return a pipeBase.Struct with fields:
        - exposure: the exposure after application of ISR
        """
        self.log.log(self.log.INFO, "Performing ISR on sensor %s" % (sensorRef.dataId))
        ccdExp = sensorRef.get('raw')
        ccd = cameraGeom.cast_Ccd(ccdExp.getDetector())

        ccdExp = isr.floatImageFromInt(ccdExp)
        
        for amp in ccd:
            ampExp = ccdExp.Factory(ccdExp, amp.getDiskAllPixels())
            amp = cameraGeom.cast_Amp(ampExp.getDetector())
    
            self.saturationDetection(ampExp, amp)

            self.overscanCorrection(ampExp, amp)

        self.biasCorrection(ccdExp, sensorRef)
        
        self.darkCorrection(ccdExp, sensorRef)
        
        self.updateVariance(ccdExp, ccd)
        
        self.flatCorrection(ccdExp, sensorRef)
        
        ccdExp = self.assembleCcd.run(ccdExp).exposure
        ccd = cameraGeom.cast_Ccd(ccdExp.getDetector())
        
        self.maskAndInterpDefect(ccdExp, ccd)
        
        self.saturationInterpolation(ccdExp)
        
        self.maskAndInterpNan(ccdExp)

        if self.config.doWrite:
            sensorRef.put(ccdExp, "visitCCD")
        
        self.display("visitCCD", ccdExp)

        return pipeBase.Struct(
            exposure = ccdExp,
        )

    def checkIsAmp(self, detector):
        """Check if a detector is of type cameraGeom.Amp

        @param Detector cameraGeom.Detector to be checked
        @return True if Amp, else False
        """
        return isinstance(detector, cameraGeom.Amp)
    
    def biasCorrection(self, exposure, dataRef):
        biasMI = sensorRef.get("bias").getMaskedImage()
        isr.biasCorrection(ccdExp.getMaskedImage(), biasMI)

    def darkCorrection(self, exposure, dataRef):
        darkexposure = dataRef['dark']
        darkScale = darkexposure.getCalib().getExptime()
        expScale = exposure.getCalib().getExptime()
        
        isr.darkCorrection(exposure.getMaskedImage(), darkexposure.getMaskedImage(), expScale, darkScale)
    
    def updateVariance(self, exposure, ccd):
        """Set the variance plane
        """
        isr.updateVariance(ccdExp.getMaskedImage(), ccd.getElectronicParams().getGain(),
            ccd.getElectronicParams().getReadNoise())

    def flatCorrection(self, exposure, dataRef):
        flatfield = dataRef['flat']

        isr.flatCorrection(
            maskedImage = exposure.getMaskedImage(),
            flatMaskedImage = flatfield.getMaskedImage(),
            scalingType = self.config.flatScalingType,
            userScale = self.config.flatUserScale,
        )

    def saturationCorrection(self, exposure, detector):
        """Perform saturation detection and interpolation as a single step.

        @warning Only suitable for CCDs with a single amplifier; otherwise you don't correctly handle
        bleed trails that go across amplifier boundaries. For more complex CCDs perform saturation detection
        first (before linearization), then CCD assembly, then saturation interpolation.
        """
        if not self.checkIsAmp(detector):
            raise RuntimeError("This method must be executed on an amp.")
        fwhm = self.config.fwhm
        grow = self.config.growSaturationFootprintSize
        maskName = self.config.saturatedMaskName
        ep = detector.getElectronicParams()
        satvalue = ep.getSaturationLevel()
        isr.saturationCorrection(
            maskedImage = exposure.getMaskedImage(),
            saturation = satvalue,
            fwhm = fwhm,
            growFootprints = grow,
            maskName = maskName,
        )
        return exposure

    def saturationDetection(self, exposure, detector):
        if not self.checkIsAmp(detector):
            raise RuntimeError("This method must be executed on an amp.")
        ep = detector.getElectronicParams()
        satvalue = ep.getSaturationLevel()
        maskName = self.config.saturatedMaskName
        isr.makeThresholdMask(
            maskedImage = exposure.getMaskedImage(),
            saturation = satvalue,
            growFootprints = 0,
            maskName = maskName,
        )
        return exposure

    def saturationInterpolation(self, exposure):
        #Don't loop over amps since saturation can cross amp boundaries
        maskName = self.config.saturatedMaskName
        fwhm = self.config.fwhm
        grow = self.config.growSaturationFootprintSize
        if self.transposeForInterpolation:
            mi = isr.transposeMaskedImage(exposure.getMaskedImage())
            isr.interpolateFromMask(mi, fwhm, growFootprints=grow, maskName=maskName)
            mi = isr.transposeMaskedImage(mi)
            exposure.setMaskedImage(mi)
        else:
            isr.interpolateFromMask(exposure.getMaskedImage(), fwhm, growFootprints=grow, maskName=maskName)
    
    def maskAndInterpDefect(self, exposure, detector):
        #Don't loop over amps since defects could cross amp boundaries
        fwhm = self.config.fwhm
        grow = self.config.growDefectFootprintSize
        defectBaseList = detector.getDefects()
        defectList = measAlg.DefectListT()
        #mask bad pixels in the camera class
        #create master list of defects and add those from the camera class
        for d in defectBaseList:
            bbox = d.getBBox()
            nd = measAlg.Defect(bbox)
            defectList.append(nd)
        isr.maskPixelsFromDefectList(exposure.getMaskedImage(), defectList, maskName='BAD')
        defectList = isr.getDefectListFromMask(exposure.getMaskedImage(), maskName='BAD', growFootprints=grow)
        if self.transposeForInterpolation:
            mi = isr.transposeMaskedImage(exposure.getMaskedImage())
            defectList = isr.transposeDefectList(defectList)
            isr.interpolateDefectList(mi, defectList, fwhm)
            mi = isr.transposeMaskedImage(mi)
            exposure.setMaskedImage(mi)
        else:
            isr.interpolateDefectList(exposure.getMaskedImage(), defectList, fwhm)

    def maskAndInterpNan(self, exposure):
        #Don't loop over amps since nans could cross amp boundaries
        fwhm = self.config.fwhm
        grow = self.config.growDefectFootprintSize
        #find unmasked bad pixels and mask them
        exposure.getMaskedImage().getMask().addMaskPlane("UNMASKEDNAN") 
        unc = UnmaskedNanCounterF()
        unc.apply(exposure.getMaskedImage())
        nnans = unc.getNpix()
        self.metadata.set("NUMNANS", nnans)
        if not nnans == 0:
		raise RuntimeError("There were %i unmasked NaNs"%(nnans))
        #get footprints of bad pixels not in the camera class
        if nnans > 0:
            undefects = isr.getDefectListFromMask(exposure.getMaskedImage(), maskName='UNMASKEDNAN', growFootprints=grow)
            #interpolate all bad pixels
            if self.transposeForInterpolation:
                mi = isr.transposeMaskedImage(exposure.getMaskedImage())
                defectList = isr.transposeDefectList(uudefects)
                isr.interpolateDefectList(mi, defectList, fwhm)
                mi = isr.transposeMaskedImage(mi)
                exposure.setMaskedImage(mi)
            else:
                isr.interpolateDefectList(exposure.getMaskedImage(), uudefects, fwhm)

    def overscanCorrection(self, exposure, detector):
        if not self.checkIsAmp(detector):
            raise RuntimeError("This method must be executed on an amp.")
        fittype = self.config.overscanFitType
        polyorder = self.config.overscanPolyOrder
        expImage = exposure.getMaskedImage().getImage()
        overscan = expImage.Factory(expImage, detector.getDiskBiasSec())
        isr.overscanCorrection(exposure.getMaskedImage(), overscan, fittype=fittype, polyorder=polyorder,
                                        imageFactory=overscan.Factory)
