import lsst.afw.image as afwImage
import lsst.pipe.base as pipeBase
from lsst.afw.detection import FootprintSet, Threshold
from lsst.pex.config import Config

__all__ = ["ProportionalLinearizeTask"]

## \addtogroup LSST_task_documentation
## \{
## \page ProportionalLinearizeTask
## \ref ProportionalLinearizeTask_ "ProportionalLinearizeTask"
## \copybrief ProportionalLinearizeTask
## \}

class ProportionalLinearizeTask(pipeBase.Task):
    """!Correct for detector non-linearity using the HSC "proportional" model

    @anchor ProportionalLinearizeTask_
    
    @section ip_isr_proportionalNonlinearity_Contents  Contents

     - @ref ip_isr_proportionalNonlinearity_Purpose
     - @ref ip_isr_proportionalNonlinearity_Initialize
     - @ref ip_isr_proportionalNonlinearity_IO
     - @ref ip_isr_proportionalNonlinearity_Config
     - @ref ip_isr_proportionalNonlinearity_Debug

    @section ip_isr_proportionalNonlinearity_Purpose  Description

    Correct for detector non-linearity using a proportional model whose coefficients are:
    - coeff[0] sqareCoeff: corrIm = uncorrIm + sqareCoeff*uncorrIm^2
    - coeff[1] must be 0 (this coefficient is called "threshold", but is unused beyond it is 0)
    - coeff[2] maxUncorr: if maxUncorr > 0 then flag pixels as SUSPECT where uncorrIm > maxUncorr
    additional coefficients are ignored (unfortunately some cameras have them)

    This is a task primarily to easily allow using other algorithms to correct non-linearity.

    @section ip_isr_proportionalNonlinearity_Initialize  Task initialisation

    @copydoc \_\_init\_\_

    @section ip_isr_proportionalNonlinearity_IO  Invoking the Task

    The main method is `run`, which takes a single exposure and applies non-linearity correction
    to it, modifying the exposure in place.

    @section ip_isr_proportionalNonlinearity_Config  Configuration parameters

    ProportionalLinearizeTask has no configuration parameters

    @section ip_isr_proportionalNonlinearity_Debug  Debug variables

    ProportionalLinearizeTask has no debug output
    """
    ConfigClass = Config
    _DefaultName = "correctNonLinearity"

    @pipeBase.timeMethod
    def run(self, exposure):
        """Correct for non-linearity using a polynomial model

        @param exposure Exposure to process

        @throw RuntimeError if there are not exactly 3 non-linearity coefficients or the second coefficient
            is non-zero.
        """
        detector = exposure.getDetector()

        if not self.doCorrect(detector):
            self.log.logdebug("Non-linearity correction not wanted for detector %s" % (detector.getId(),))
            return

        didCorrect = False  # did we apply linearity correction to any amplifier?
        for amp in detector:
            squareCoeff, coeff1, maxUncorr = amp.getLinearityCoeffs()[0:3]

            if coeff1 != 0.0:
                raise RuntimeError(
                    "coeff[1] must be 0 for PROPORTIONAL non-linearity correction; saw %g "
                     "for detector %r amplifier %r" % (coeff1, detector.getId(), amp.getName()))
                if (squareCoeff, maxUncorr) != (squareCoeff, maxUncorr):
                    raise RuntimeError("One or more coefficients not finite: %s" %
                        (squareCoeff, coeff1, maxUncorr))

            if maxUncorr <= 0 and squareCoeff == 0.0:
                continue  # nothing to do

            ampImage = afwImage.MaskedImageF(exposure.getMaskedImage(), amp.getBBox(), afwImage.PARENT)

            if maxUncorr > 0:
                FootprintSet(ampImage, Threshold(maxUncorr), "SUSPECT")

            if squareCoeff != 0:
                ampArr = ampImage.getImage().getArray()
                ampArr *= 1.0 + squareCoeff*ampArr

            didCorrect = True

        if didCorrect:
            self.log.log(self.log.INFO, "Applied linearity corrections to detector %s" % (detector.getId()))

    def doCorrect(self, detector):
        """!Return True if non-linearity correction is wanted for this detector

        Return True if linearity type is "PROPORTIONAL", False if "NONE", raise an exception otherwise.

        @param[in] detector  detector info (an lsst.afw.cameraGeom.Detector)

        @note only checks the first amplifier; if the detector requires different kinds of
        non-linearity correction for different amplifiers then use a different task.

        @throw RuntimeError if correction type is not "NONE" or "PROPORTIONAL"
        """
        amp = detector[0]

        linearityType = amp.getLinearityType()
        if linearityType == "NONE":
            return False
        elif linearityType == "PROPORTIONAL":
            numCoeffs = len(amp.getLinearityCoeffs())
            if numCoeffs < 3:
                raise RuntimeError("Need 3 linearity coefficients; only found %d" % (numCoeffs,))
            return True

        raise RuntimeError("Unsupported linearity type %r for detector %r amplifier %r; "
            "must be PROPORTIONAL or NONE" % (linearityType, detector.getId(), amp.getName()))
