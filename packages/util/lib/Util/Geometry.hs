{-# LANGUAGE FlexibleInstances, FlexibleContexts #-}
--, TypeSynonymInstances,  #-}
{-# LANGUAGE TypeFamilies #-}
-----------------------------------------------------------------------------
{- |
Module      :  Util.Geometry
Copyright   :  (c) Alberto Ruiz 2006-12
License     :  GPL

Maintainer  :  Alberto Ruiz (aruiz at um dot es)
Stability   :  provisional

Projective geometry utilities.

-}
-----------------------------------------------------------------------------

module Util.Geometry
(
  -- * Basic Types
    Point(..), HPoint(..), HLine(..),
    Point3D(..), HPoint3D(..), HLine3D(..), HPlane(..), 

    Homography, Camera, Homography3D,

    Conic, DualConic, Quadric, DualQuadric,

    Vectorable(toVector), Matrixlike(toMatrix),
    mkTrans, Row2(..),Row3(..),Row4(..),Col2(..),Col3(..),Col4(..),

  -- * Transformations

  -- * Util
    joinPoints, meetLines

) where

import Util.Misc(Mat,Vec)
import Numeric.LinearAlgebra(fromList,(@>),(><),fromRows,(<>))

----------------------------------------------------------------------

class Vectorable a
  where
    toVector   :: a -> Vec
    fromVector :: Vec -> a

class Matrixlike a where
    toMatrix :: a -> Mat
    unsafeFromMatrix :: Mat -> a

----------------------------------------------------------------------

-- | inhomogenous 2D point
data Point = Point {px :: !Double, py :: !Double} deriving (Eq, Show, Read)

instance Vectorable Point where
    toVector (Point x y) = fromList [x,y]
    fromVector v = Point (v@>0) (v@>1)


-- | inhomogenous 2D point
data HPoint = HPoint !Double !Double !Double deriving (Eq, Show, Read)

instance Vectorable HPoint where
    toVector (HPoint x y w) = fromList [x,y,w]
    fromVector v = HPoint (v@>0) (v@>1) (v@>2)


-- | inhomogenous 3D point
data Point3D = Point3D !Double !Double !Double deriving (Eq, Show, Read)

instance Vectorable Point3D where
    toVector (Point3D x y z) = fromList [x,y,z]
    fromVector v = Point3D (v@>0) (v@>1) (v@>3)


-- | homogenous 3D point
data HPoint3D = HPoint3D !Double !Double !Double !Double deriving (Eq, Show, Read)

instance Vectorable HPoint3D where
    toVector (HPoint3D x y z w) = fromList [x,y,z,w]
    fromVector v = HPoint3D (v@>0) (v@>1) (v@>2) (v@>3)


-- | 2D line
data HLine = HLine {aLn, bLn, cLn :: !Double} deriving (Eq, Show, Read)

instance Vectorable HLine where
    toVector (HLine a b c) = fromList [a,b,c]
    fromVector v = HLine (v@>0) (v@>1) (v@>3)


-- | 3D line (provisional)
newtype HLine3D = HLine3D Mat deriving (Eq, Show, Read)

instance Matrixlike HLine3D where
    toMatrix (HLine3D m) = m
    unsafeFromMatrix = HLine3D

-- | 3D plane
data HPlane = HPlane !Double !Double !Double !Double deriving (Eq, Show, Read)

instance Vectorable HPlane where
    toVector (HPlane a b c d) = fromList [a,b,c,d]
    fromVector v = HPlane (v@>0) (v@>1) (v@>2) (v@>3)


-- | projective transformation P2->P2
newtype Homography = Homography Mat deriving (Eq, Show, Read)

instance Matrixlike Homography where
    toMatrix (Homography m) = m
    unsafeFromMatrix = Homography

-- | projective transformation P3->P2
newtype Camera = Camera Mat deriving (Eq, Show, Read)

instance Matrixlike Camera where
    toMatrix (Camera m) = m
    unsafeFromMatrix = Camera

-- | projective transformation P3->P3
newtype Homography3D = Homography3D Mat deriving (Eq, Show, Read)

instance Matrixlike Homography3D where
    toMatrix (Homography3D m) = m
    unsafeFromMatrix = Homography3D

newtype Conic = Conic Mat deriving (Eq, Show, Read)

instance Matrixlike Conic where
    toMatrix (Conic m) = m
    unsafeFromMatrix = Conic


newtype Quadric = Quadric Mat deriving (Eq, Show, Read)

instance Matrixlike Quadric where
    toMatrix (Quadric m) = m
    unsafeFromMatrix = Quadric

newtype DualConic = DualConic Mat deriving (Eq, Show, Read)

instance Matrixlike DualConic where
    toMatrix (DualConic m) = m
    unsafeFromMatrix = DualConic

newtype DualQuadric = DualQuadric Mat deriving (Eq, Show, Read)

instance Matrixlike DualQuadric where
    toMatrix (DualQuadric m) = m
    unsafeFromMatrix = DualQuadric


data Dim2 a = Dim2 !a !a
data Dim3 a = Dim3 !a !a !a
data Dim4 a = Dim4 !a !a !a !a



data Col2 = Col2 !Double !Double
data Col3 = Col3 !Double !Double !Double
data Col4 = Col4 !Double !Double !Double !Double

data Row2 a = Row2 a a
data Row3 a = Row3 a a a
data Row4 a = Row4 a a a a

instance Vectorable Col2 where
    toVector (Col2 x1 x2) = fromList [x1,x2]
    fromVector v = Col2 (v@>0) (v@>1)

instance Vectorable Col3 where
    toVector (Col3 x1 x2 x3) = fromList [x1,x2,x3]
    fromVector v = Col3 (v@>0) (v@>1) (v@>2)

instance Vectorable Col4 where
    toVector (Col4 x1 x2 x3 x4) = fromList [x1,x2,x3,x4]
    fromVector v = Col4 (v@>0) (v@>1) (v@>2) (v@>3)


matrix2x2 :: Row2 Col2 -> Mat
matrix2x2 (Row2 (Col2 x1 x2)
                (Col2 x3 x4) ) = (3><3) [x1,x2,
                                         x3,x4]

matrix3x3 :: Row3 Col3 -> Mat
matrix3x3 (Row3 (Col3 x1 x2 x3)
                (Col3 x4 x5 x6)
                (Col3 x7 x8 x9) ) = (3><3) [x1,x2,x3,
                                            x4,x5,x6,
                                            x7,x8,x9]

matrix3x4 :: Row3 Col4 -> Mat
matrix3x4 (Row3 r1 r2 r3) = fromRows (map toVector [r1,r2,r3])

matrix4x4 :: Row4 Col4 -> Mat
matrix4x4 (Row4 r1 r2 r3 r4) = fromRows (map toVector [r1,r2,r3,r4])

type family MatrixShape  (m :: *)

type instance MatrixShape Homography = Row3 Col3
type instance MatrixShape Camera = Row3 Col4
type instance MatrixShape Homography3D = Row4 Col4
type instance MatrixShape Conic = Row3 Col3
type instance MatrixShape DualConic = Row3 Col3
type instance MatrixShape Quadric = Row4 Col4
type instance MatrixShape DualQuadric = Row4 Col4

class MatrixElem t where
    fromElements :: t -> Mat

instance MatrixElem (Row3 Col3) where
    fromElements = matrix3x3

instance MatrixElem (Row2 Col2) where
    fromElements = matrix2x2

instance MatrixElem (Row3 Col4) where
    fromElements = matrix3x4

instance MatrixElem (Row4 Col4) where
    fromElements = matrix4x4

mkTrans :: (Matrixlike t, MatrixElem (MatrixShape t)) => MatrixShape t -> t
mkTrans = unsafeFromMatrix . fromElements


crossMat v = (3><3) [ 0,-c, b,
                   c, 0,-a,
                  -b, a, 0]
    where a = v@>0
          b = v@>1
          c = v@>2

joinPoints :: HPoint -> HPoint -> HLine
joinPoints p q = HLine (v@>0) (v@>1) (v@>2) where v = crossMat (toVector p) <> (toVector q)

meetLines :: HLine -> HLine -> HPoint
meetLines l m = HPoint (v@>0) (v@>1) (v@>2) where v = crossMat (toVector l) <> (toVector m)


