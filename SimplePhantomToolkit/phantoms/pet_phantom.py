#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 14 11:19:09 2017

@author: marcel
"""

from SimplePhantomToolkit import Cylinder, geometry


class PetPhantom(Cylinder):
    LENGTH = 200
    DIAMETER = 200
    VOLUME = 6283
    def __init__(self, image, margin = 0.7):
        """ Homogeneous PET phantom for Siemens Camera's. """

        # assume the centroid of the image is the center of the physical phantom
        anchor = geometry.centroid(image)

        axis = [[0, 0, -1 * margin * 0.5 * self.LENGTH],
                [0, 0, margin * 0.5 * self.LENGTH]]

        radius = self.DIAMETER/2 * margin

        super(PetPhantom, self).__init__(image = image, axis = axis,
             radius = radius, anchor = anchor)

        self.reset()

