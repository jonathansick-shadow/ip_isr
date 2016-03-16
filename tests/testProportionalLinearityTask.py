#!/usr/bin/env python
from __future__ import absolute_import, division, print_function
#
# LSST Data Management System
# Copyright 2008-2013 LSST Corporation.
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
import unittest

import numpy as np

import lsst.utils
import lsst.utils.tests
import lsst.afw.geom as afwGeom
from lsst.afw.image import ExposureF
from lsst.afw.cameraGeom.testUtils import DetectorWrapper
from lsst.ip.isr import ProportionalLinearizeTask

def makeSimpleDetector():
    """Make a very minimal detector with four amplifiers arranged 2x2 whose bboxes covers the full image
    """
    ampDim = afwGeom.Extent2I(500, 300)
    expBBox = afwGeom.Box2I(afwGeom.Point2I(0, 0), ampDim*2)
    dw = DetectorWrapper(bbox=expBBox, numAmps=4)
    detector = dw.detector
    for i, amp in enumerate(detector):
        ampMin = afwGeom.Point2I(
            ampDim[0] if i in (1, 3) else 0,
            ampDim[1] if i in (2, 3) else 0,
        )
        amp.setBBox(afwGeom.Box2I(ampMin, ampDim))
    return detector


class ProportionalLinearizeTestCase(lsst.utils.tests.TestCase):
    def setUp(self):
        # make a detector object and a blank exposure containing it
        detector = makeSimpleDetector()
        self.exposure = ExposureF(detector.getBBox())
        self.exposure.setDetector(detector)

    def tearDown(self):
        del self.exposure

    def testBasicCorrection(self):
        """test ProportionalLinearizeTask
        """
        suspectValue = 50000 # a largish value smaller than saturation
        squareCoeff = 1e-5

        # compute amplifier data region as follows:
        # - image and variance to a ramp that goes from 0 to a value larger than suspectValue
        # - mask is zero (already true)
        detector = self.exposure.getDetector()

        imArr, maskArr, varArr = self.exposure.getMaskedImage().getArrays()
        imArr.flat[:] = np.linspace(0, suspectValue + 1000, imArr.size)
        begImArr = np.copy(imArr)
        maskArr[:] = 0 # mask
        varArr[:] = imArr[:]

        desImArr = imArr*(1.0 + squareCoeff*imArr)
        suspectMaskValue = self.exposure.getMaskedImage().getMask().getPlaneBitMask(["SUSPECT"])
        desMaskArr = np.where(imArr > suspectValue, suspectMaskValue, 0)
        desArrSet = (desImArr, desMaskArr, begImArr)

        # set appropriate linearity configuration
        # and set data pixels as follows:
        # 
        detector = self.exposure.getDetector()
        for amp in detector:
            amp.setLinearityType("PROPORTIONAL")
            amp.setLinearityCoeffs((squareCoeff, 0, suspectValue))

        task = ProportionalLinearizeTask()
        task.run(self.exposure)
        endMi = self.exposure.getMaskedImage()

        self.assertMaskedImagesNearlyEqual(endMi, desArrSet)

    def testNoCorrection(self):
        begMi = self.exposure.getMaskedImage()
        begMi = begMi.Factory(begMi, True) # deep copy
        detector = self.exposure.getDetector()
        for amp in detector:
            amp.setLinearityType("NONE")
        task = ProportionalLinearizeTask()
        task.run(self.exposure)
        endMi = self.exposure.getMaskedImage()
        self.assertMaskedImagesNearlyEqual(begMi, endMi, rtol=0, atol=0)

def suite():
    lsst.utils.tests.init()
    suites = []
    suites += unittest.makeSuite(ProportionalLinearizeTestCase)
    suites += unittest.makeSuite(lsst.utils.tests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(shouldExit=False):
    lsst.utils.tests.run(suite(), shouldExit)

if __name__ == "__main__":
    run(True)
