{-# LANGUAGE FlexibleInstances, FlexibleContexts #-}
--, TypeSynonymInstances,
{-# LANGUAGE TypeOperators #-}
{-# LANGUAGE MultiParamTypeClasses #-}
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

    Vectorlike(..), Matrixlike(..), Tensorial(..),
    mkTrans, Dim2(..),Dim3(..),Dim4(..),

  -- * Transformations

    Transformable(..), Composable(..),

  -- * Util
    joinPoints, meetLines,
    
  -- * Conversion
    Inhomog(..), Homog(..)

) where

import Util.Misc(Mat,Vec)
import Numeric.LinearAlgebra(fromList,(@>),(><),toRows,fromRows,(<>),trans,inv)
--import Foreign.Storable(Storable)
import qualified Numeric.LinearAlgebra.Tensor as T
import qualified Numeric.LinearAlgebra.Array.Util as UT

----------------------------------------------------------------------

class Vectorlike a
  where
    toVector   :: a -> Vec
    fromVector :: Vec -> a

class Matrixlike a where
    toMatrix :: a -> Mat
    unsafeFromMatrix :: Mat -> a

----------------------------------------------------------------------

-- | inhomogenous 2D point
data Point = Point {px :: !Double, py :: !Double} deriving (Eq, Show, Read)

instance Vectorlike Point where
    toVector (Point x y) = fromList [x,y]
    fromVector v = Point (v@>0) (v@>1)


-- | inhomogenous 2D point
data HPoint = HPoint !Double !Double !Double deriving (Eq, Show, Read)

instance Vectorlike HPoint where
    toVector (HPoint x y w) = fromList [x,y,w]
    fromVector v = HPoint (v@>0) (v@>1) (v@>2)


-- | inhomogenous 3D point
data Point3D = Point3D !Double !Double !Double deriving (Eq, Show, Read)

instance Vectorlike Point3D where
    toVector (Point3D x y z) = fromList [x,y,z]
    fromVector v = Point3D (v@>0) (v@>1) (v@>2)


-- | homogenous 3D point
data HPoint3D = HPoint3D !Double !Double !Double !Double deriving (Eq, Show, Read)

instance Vectorlike HPoint3D where
    toVector (HPoint3D x y z w) = fromList [x,y,z,w]
    fromVector v = HPoint3D (v@>0) (v@>1) (v@>2) (v@>3)


-- | 2D line
data HLine = HLine {aLn, bLn, cLn :: !Double} deriving (Eq, Show, Read)

instance Vectorlike HLine where
    toVector (HLine a b c) = fromList [a,b,c]
    fromVector v = HLine (v@>0) (v@>1) (v@>2)


-- | 3D line (provisional)
newtype HLine3D = HLine3D Mat deriving (Eq, Show, Read)

instance Matrixlike HLine3D where
    toMatrix (HLine3D m) = m
    unsafeFromMatrix = HLine3D

-- | 3D plane
data HPlane = HPlane !Double !Double !Double !Double deriving (Eq, Show, Read)

instance Vectorlike HPlane where
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

instance Vectorlike (Dim2 Double) where
    toVector (Dim2 x1 x2) = fromList [x1,x2]
    fromVector v = Dim2 (v@>0) (v@>1)

instance Vectorlike (Dim3 Double) where
    toVector (Dim3 x1 x2 x3) = fromList [x1,x2,x3]
    fromVector v = Dim3 (v@>0) (v@>1) (v@>2)

instance Vectorlike (Dim4 Double) where
    toVector (Dim4 x1 x2 x3 x4) = fromList [x1,x2,x3,x4]
    fromVector v = Dim4 (v@>0) (v@>1) (v@>2) (v@>3)


matrix2x2 :: Dim2 (Dim2 Double) -> Mat
matrix2x2 (Dim2 (Dim2 x1 x2)
                (Dim2 x3 x4) ) = (3><3) [x1,x2,
                                         x3,x4]

matrix3x3 :: Dim3 (Dim3 Double) -> Mat
matrix3x3 (Dim3 (Dim3 x1 x2 x3)
                (Dim3 x4 x5 x6)
                (Dim3 x7 x8 x9) ) = (3><3) [x1,x2,x3,
                                            x4,x5,x6,
                                            x7,x8,x9]

matrix3x4 :: Dim3 (Dim4 Double) -> Mat
matrix3x4 (Dim3 r1 r2 r3) = fromRows (map toVector [r1,r2,r3])

matrix4x4 :: Dim4 (Dim4 Double) -> Mat
matrix4x4 (Dim4 r1 r2 r3 r4) = fromRows (map toVector [r1,r2,r3,r4])

type family MatrixShape  (m :: *)

type Dim2x2 = Dim2 (Dim2 Double)
type Dim3x3 = Dim3 (Dim3 Double)
type Dim3x4 = Dim3 (Dim4 Double)
type Dim4x4 = Dim4 (Dim4 Double)

type instance MatrixShape Homography = Dim3x3
type instance MatrixShape Camera = Dim3x4
type instance MatrixShape Homography3D = Dim4x4
type instance MatrixShape Conic = Dim3x3
type instance MatrixShape DualConic = Dim3x3
type instance MatrixShape Quadric = Dim4x4
type instance MatrixShape DualQuadric = Dim4x4

class MatrixElem t where
    fromElements :: t -> Mat

instance MatrixElem Dim3x3 where
    fromElements = matrix3x3

instance MatrixElem Dim2x2 where
    fromElements = matrix2x2

instance MatrixElem Dim3x4 where
    fromElements = matrix3x4

instance MatrixElem Dim4x4 where
    fromElements = matrix4x4

mkTrans :: (Matrixlike t, MatrixElem (MatrixShape t)) => MatrixShape t -> t
mkTrans = unsafeFromMatrix . fromElements

crossMat :: Vec -> Mat
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

class Transformable t x
  where
    type TResult t x :: *
    apTrans :: t -> x -> TResult t x
    infixl 2 <|
    (<|) :: t -> x -> TResult t x
    (<|) = apTrans
    infixl 2 ◁  -- 25c1
    (◁) :: t -> x -> TResult t x
    (◁) = apTrans


instance Transformable Homography [HPoint]
  where
    type TResult Homography [HPoint] = [HPoint]
    apTrans h = (map fromVector . toRows) . (<> trans (toMatrix h)) . fromRows . (map toVector)

instance Transformable Homography [Point]
  where
    type TResult Homography [Point] = [Point]
    apTrans h = map inhomog . apTrans h . map homog -- FIXME

instance Transformable Homography [HLine]
  where
    type TResult Homography [HLine] = [HLine]
    apTrans h = (map fromVector . toRows) . (<> inv (toMatrix h)) . fromRows . (map toVector)


class Composable s t
  where
    type s :.: t :: *
    compTrans :: s -> t -> s :.: t  -- s . t
    (·) :: s -> t -> s :.: t
    infixr 8  · -- utf8 2299
    (·) = compTrans
    (⊙) :: s -> t -> s :.: t
    infixr 8  ⊙ -- utf8 2299
    (⊙) = compTrans

instance Composable Homography Homography
  where
    type Homography :.: Homography = Homography
    compTrans s t = unsafeFromMatrix (toMatrix s <> toMatrix t)

instance Composable Homography3D Homography3D
  where
    type Homography3D :.: Homography3D = Homography3D
    compTrans s t = unsafeFromMatrix (toMatrix s <> toMatrix t)

instance Composable Camera Homography3D
  where
    type Camera :.: Homography3D = Camera
    compTrans s t = unsafeFromMatrix (toMatrix s <> toMatrix t)

instance Composable Homography Camera
  where
    type Homography :.: Camera = Camera
    compTrans s t = unsafeFromMatrix (toMatrix s <> toMatrix t)



class Inhomog x
  where
    type HResult x :: *
    homog :: x -> HResult x

instance Inhomog Point
  where
    type HResult Point = HPoint
    homog (Point x y) = HPoint x y 1

instance Inhomog Point3D
  where
    type HResult Point3D = HPoint3D
    homog (Point3D x y z) = HPoint3D x y z 1



class Homog x
  where
    type IHResult x :: *
    inhomog :: x -> IHResult x

instance Homog HPoint
  where
    type IHResult HPoint = Point
    inhomog (HPoint x y w) = Point (x/w) (y/w)

instance Homog HPoint3D
  where
    type IHResult HPoint3D = Point3D
    inhomog (HPoint3D x y z w) = Point3D (x/w) (y/w) (z/w)


class Tensorial x
  where
    toTensor :: x -> T.Tensor Double


instance Tensorial HPoint
  where
    toTensor (HPoint x y w) = T.vector [x,y,w]
    
instance Tensorial HPoint3D
  where
    toTensor (HPoint3D x y z w) = T.vector [x,y,z,w]

instance Tensorial HLine
  where
    toTensor (HLine a b c) = T.covector [a,b,c]

instance Tensorial HPlane
  where
    toTensor (HPlane a b c d) = T.covector [a,b,c,d]

instance Tensorial HLine3D
  where
    toTensor (HLine3D m) = UT.fromMatrix T.Contra T.Contra m

