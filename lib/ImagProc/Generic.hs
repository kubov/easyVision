{-# OPTIONS_GHC -fglasgow-exts #-}

-----------------------------------------------------------------------------
{- |
Module      :  ImagProc.Generic
Copyright   :  (c) Alberto Ruiz 2007
License     :  GPL-style

Maintainer  :  Alberto Ruiz (aruiz at um dot es)
Stability   :  very provisional
Portability :  hmm...

-}
-----------------------------------------------------------------------------

module ImagProc.Generic (
  GImg(..)
, blockImage
, warp, warpOn
, constImage, cloneClear
, Channels(..), channels,
  channelsFromRGB
)
where

import ImagProc.Ipp
import Foreign.C.Types(CUChar)
import Foreign hiding (shift)
import Numeric.LinearAlgebra
import Vision

class Image image => GImg pixel image | pixel -> image, image -> pixel where
    zeroP :: pixel
    set :: pixel -> ROI -> image -> IO ()
    copy :: image -> ROI -> image -> ROI -> IO ()
    -- | Resizes the roi of a given image.
    resize :: Size -> image -> image
    clone :: image -> IO image
    warpOnG :: [[Double]] -> image -> image -> IO ()
    fromYUV :: ImageYUV -> image
    toYUV :: image -> ImageYUV

----------------------------------

instance GImg CUChar ImageGray where
    zeroP = 0
    set = set8u
    copy = copyROI8u
    resize = resize8u
    clone = ioCopy_8u_C1R id
    warpOnG = warpOn8u
    fromYUV = yuvToGray
    toYUV = grayToYUV

instance GImg Float ImageFloat where
    zeroP = 0
    set = set32f
    copy = copyROI32f
    resize = resize32f
    clone = ioCopy_32f_C1R id
    warpOnG = warpOn32f
    fromYUV = float . yuvToGray
    toYUV = grayToYUV . toGray

instance GImg (CUChar,CUChar,CUChar) ImageRGB where
    zeroP = (0,0,0)
    set (r,g,b) = set8u3 r g b
    copy = copyROI8u3
    resize = resize8u3
    clone = ioCopy_8u_C3R id
    warpOnG = warpOn8u3
    fromYUV = yuvToRGB
    toYUV = rgbToYUV



-- modifies the undefined region of an image.
clearNoROI :: Image a => (ROI -> a -> IO ()) -> a -> IO ()
clearNoROI fun im = mapM_ ((flip fun) im) (invalidROIs im)

--cloneClear' :: Image a => a -> IO ()
cloneClear im = do
    r <- clone im
    clearNoROI (set zeroP) r
    return r

--------------------------------------------------------------------------

-- | joins images
blockImage :: GImg p img  => [[img]] -> img
blockImage = columnImage . map rowImage

--rowImage :: [ImageGray] -> ImageGray
rowImage l = unsafePerformIO $ do
    let r = maximum (map (height.size) l)
        c = maximum (map (width.size) l)
        n = length l
    res <- image (Size r (c*n))
    let roi0 = theROI (head l)
        rois = take n $ iterate (shift (0,c)) roi0
        f r i = copy i roi0 res r
    sequence_ $ zipWith f rois l
    return res

--columnImage :: [ImageGray] -> ImageGray
columnImage l = unsafePerformIO $ do
    let r = maximum (map (height.size) l)
        c = maximum (map (width.size) l)
        n = length l
    res <- image (Size (r*n) c)
    let roi0 = theROI (head l)
        rois = take n $ iterate (shift (r,0)) roi0
        f r i = copy i roi0 res r
    sequence_ $ zipWith f rois l
    return res


--------------------------------------------------------------------

adapt dst h src = toLists $ inv (pixelToPointTrans (size dst)) <> h <> pixelToPointTrans (size src)

-- | Apply a homography (defined on normalized points, see 'pixelsToPoints') to an image.
{-
warp :: Size              -- ^ desired size of the result
     -> Matrix Double     -- ^ homography
     -> ImageFloat        -- ^ source image
     -> IO ImageFloat     -- ^ result
-}

warp out s h im = unsafePerformIO $ do
    r <- image s
    set out (theROI r) r
    warpOn h r im
    return r

-- | The same as 'warp', but the result is written over a preexisting image.
{-
warpOn :: Matrix Double   -- ^ homography
       -> ImageFloat      -- ^ destination image
       -> ImageFloat      -- ^ source image
       -> IO ()
-}
warpOn h r im = warpOnG (adapt r h im) r im

---------------------------------------------------------------------------------

-- hmmmmmmmm
constImage val sz = unsafePerformIO $ do
    z <- image sz
    set val (theROI z) z
    return z

---------------------------------------------------------------------------------

data Channels = CHIm
    { yuv  :: ImageYUV
    , gray :: ImageGray
    , rgb  :: ImageRGB
    , rCh  :: ImageGray
    , gCh  :: ImageGray
    , bCh  :: ImageGray
    , hsv  :: ImageRGB
    , hCh  :: ImageGray
    , sCh  :: ImageGray
    }

channels :: ImageYUV -> Channels
channels img = CHIm
    { yuv = img
    , gray = fromYUV img
    , rgb = rgbAux
    , rCh = getChannel 0 rgbAux
    , gCh = getChannel 1 rgbAux
    , bCh = getChannel 2 rgbAux
    , hsv = hsvAux
    , hCh = getChannel 0 hsvAux
    , sCh = getChannel 1 hsvAux
    }
    where rgbAux = fromYUV img
          hsvAux = rgbToHSV rgbAux

channelsFromRGB :: ImageRGB -> Channels
channelsFromRGB img = CHIm
    { yuv = yuvAux
    , gray = fromYUV yuvAux
    , rgb = img
    , rCh = getChannel 0 img
    , gCh = getChannel 1 img
    , bCh = getChannel 2 img
    , hsv = hsvAux
    , hCh = getChannel 0 hsvAux
    , sCh = getChannel 1 hsvAux
    }
    where yuvAux = toYUV img
          hsvAux = rgbToHSV img
