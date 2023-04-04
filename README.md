Poisson Image Editing

+-----------------------+-----------------------+-----------------------+
| Patrick P´erez∗       | Michel Gangnet†       | > Andrew Blake‡       |
+=======================+=======================+=======================+
+-----------------------+-----------------------+-----------------------+

Microsoft Research UK

Abstract

Using generic interpolation machinery based on solving Poisson
equations, a variety of novel tools are introduced for seamless edit-ing
of image regions. The first set of tools permits the seamless
importation of both opaque and transparent source image regions into a
destination region. The second set is based on similar math-ematical
ideas and allows the user to modify the appearance of the image
seamlessly, within a selected region. These changes can be arranged to
affect the texture, the illumination, and the color of ob-jects lying in
the region, or to make tileable a rectangular selection.

**CR Categories:** I.3.3 \[Computer Graphics\]: Picture/Image
Generation---Display algorithms; I.3.6 \[Computer Graphics\]:
Methodology and Techniques---Interaction techniques; I.4.3 \[Im-age
Processing and Computer Vision\]: Enhancement---Filtering;

**Keywords:** Interactive image editing, image gradient, guided
in-terpolation, Poisson equation, seamless cloning, selection editing

1 Introduction

Image editing tasks concern either global changes (color/intensity
corrections, filters, deformations) or local changes confined to a
se-lection. Here we are interested in achieving local changes, ones that
are restricted to a region manually selected, in a seamless and
effortless manner. The extent of the changes ranges from slight
dis-tortions to complete replacement by novel content. Classic tools to
achieve that include image filters confined to a selection, for slight
changes, and interactive cut-and-paste with cloning tools for com-plete
replacements. With these classic tools, changes in the selected regions
result in visible seams, which can be only partly hidden, subsequently,
by feathering along the border of the selected region.

We propose here a generic machinery from which different tools for
seamless editing and cloning of a selection region can be de-rived. The
mathematical tool at the heart of the approach is the Poisson partial
differential equation with Dirichlet boundary con-ditions which
specifies the Laplacian of an unknown function over the domain of
interest, along with the unknown function values over the boundary of
the domain. The motivation is twofold.

First, it is well known to psychologists \[Land and McCann 1971\] that
slow gradients of intensity, which are suppressed by the Lapla-cian
operator, can be superimposed on an image with barely notice-

> ∗e-mail: pperez@microsoft.com\
> †e-mail:mgangnet@microsoft.com\
> ‡e-mail:ablake@microsoft.com
>
> Permission to make digital/hard copy of part of all of this work for
> personal or
>
> classroom use is granted without fee provided that the copies are not
> made or
>
> distributed for profit or commercial advantage, the copyright notice,
> the title of the
>
> publication, and its date appear, and notice is given that copying is
> by permission
>
> of ACM, Inc. To copy otherwise, to republish, to post on servers, or
> to redistribute
>
> to lists, requires prior specific permission and/or a fee.
>
> © 2003 ACM 0730-0301/03/0700-0313 \$5.00
>
> 313

![](vertopal_952755e54fbc4a5886c9aadf32ad9f47/media/image1.png){width="1.1388888888888888in"
height="0.9861111111111112in"}

object, associated edgels are removed; to add an object, associated
edgels as well as color values on both sides of each of these edgels are
incorporated. The new image is then obtained by interpolat-ing smoothly
the colors associated to the new set of edgels. This amounts to solving
a Laplace equation (a Poisson equation with a null right hand side) with
Dirichlet boundary conditions given by colors around edgels. Editing
edgels and associated colors is not always simple. In addition, image
details are lost when converting to and from the contour domain, which
might be undesirable. The sparse edgel-based representation is indeed
incomplete, as opposed to the related representation based on wavelet
extrema \[Mallat and Zhong 1992\], which are complete but less adapted
to manual edit-ing.

In \[Lewis 2001\], spots are removed from fur images by separat-ing out
the brightness component from details in a selected region and replacing
the brightness by harmonic interpolation (solving a Laplace equation) of
the brightness at the selection boundary.

In terms of image editing functionalities, two existing techniques
achieve seamless cloning as the basic instance of our system does. The
first one is Adobec⃝Photoshopc⃝7's Healing Brush \[Adobec⃝2002\]. To the
best of our knowledge, the technique used by this tool has not been
published. Therefore, we don't know whether or not it uses a Poisson
solver.

The second technique is the multiresolution image blending pro-posed in
\[Burt and Adelson 1983\]. The idea is to use a multires-olution
representation, namely a Laplacian pyramid, of the images of interest.
The content of the source image region is mixed, within each resolution
band independently, with its new surrounding in the destination image.
The final composite image is then recov-ered by adding up the different
levels of the new composite Lapla-cian pyramid thus obtained. The
technique results in multiresolu-tion mixing where finest details are
averaged very locally around the boundary of the selection, while lower
frequencies are mixed over much larger distances around these
boundaries. This fast tech-nique achieves an approximate insertion of
the source Laplacian in the destination region (on the first level of
the Laplacian pyra-mid) whereas we perform this Laplacian insertion
exactly via the solution of a Poisson equation. More importantly,
multiresolution blending incorporates data from distant source and
destination pix-els, via the upper levels of the pyramid, within the
final composite image. This long range mixing, which might be
undesirable, does not occur in our technique. In addition, our system
offers extended functionality besides opaque seamless cloning, see
Sections 3 and 4.

Finally, whereas we propose a guided interpolation framework, with the
guidance being specified by the user, e.g., in the form of a source
image in the case of seamless cloning, various interpo-lation methods
have been proposed to fill in image regions auto-matically using only
the knowledge of the boundary conditions. A first class of such
approaches is composed of inpainting tech-niques \[Ballester et al.
2001; Bertalmio et al. 2000\] where PDE-based interpolation methods are
devised such as to continue the isophotes hitting the boundary of the
selected region. The PDEs to be solved are more complex than the Poisson
equation, and work only for bridging fairly narrow gaps in relatively
texture-free re-gions. Example-based interpolation methods \[Barret and
Cheney 2002; Bornard et al. 2002; Efros and Leung 1999\] where the new
image region is synthesized using an arrangement of many small patches
are an interesting alternative to inpainting. These methods handle large
holes and textured boundaries in a more convincing way. Moreover, they
can also be used to import textures as shown in \[Efros and Freeman
2001; Hertzmann et al. 2001\].

> 314
>
> Therefore, inside Ω, the additive correction ˜*f* is a membrane
> inter-polant of the mismatch (*f*∗− *g*) between the source and the
> desti-nation along the boundary ∂Ω. This particular instance of guided
> interpolation is used for seamless cloning in Section 3.
>
> Discrete Poisson solver The variational problem (3), and the
> associated Poisson equation with Dirichlet boundary conditions (4),
> can be discretized and solved in a number of ways.

For discrete images the problem can be discretized naturally us-ing the
underlying discrete pixel grid. Without loss of generality, we will keep
the same notations for the continuous objects and their discrete
counterparts: *S*, Ω now become finite point sets defined on an infinite
discrete grid. Note that *S* can include all the pixels of an image or
only a subset of them. For each pixel *p* in *S*, let *Np* be the set of
its 4-connected neighbors which are in *S*, and let ⟨*p*,*q*⟩denote a
pixel pair such that *q* ∈ *Np*. The boundary of Ω is now∂Ω = {*p* ∈
*S*\\Ω : *Np* ∩Ω ̸= /0}. Let *fp* be the value of *f* at *p*. The task
is to compute the set of intensities *f*\|Ω =

bitrary shape, it is best to discretize the variational problem (3)
di-For Dirichlet boundary conditions defined on a boundary of ar- *fp*,
*p* ∈ Ω .

> rectly, rather than the Poisson equation (4). The finite difference
> discretization of (3) yields the following discrete, quadratic
> opti-mization problem:
>
> min *f*\|Ω⟨*p*,*q*⟩∩Ω̸=/0 ∑ (*fp* − *fq* −*vpq*)2, with *fp* =
> *f*∗*p*,for all *p* ∈ ∂Ω, (6)

where *vpq* is the projection of **v**(*[p]{.underline}*+[*q*
2]{.underline}) on the oriented edge \[*p*,*q*\], i.e., *vpq* =
**v**(*[p]{.underline}*+*[q]{.underline}* neous linear equations:
[2]{.underline})· ⃗*pq*. Its solution satisfies the following simulta-

> for all *p* ∈ Ω, \|*Np*\|*fp* − ∑*q*∈*Np*∩Ω *fq* = *q*∈*Np*∩∂Ω ∑
> *f*∗*q*+∑ *q*∈*Np* *vpq*. (7)
>
> When Ω contains pixels on the border of *S*, which happens for
> in-stance when Ω extends to the edge of the pixel grid, these pixels
> have a truncated neighborhood such that \|*Np*\| \< 4. Note that for
> pixels *p* interior to Ω, that is, *Np* ⊂ Ω, there are no boundary
> terms in the right hand side of (7), which reads:
>
> \|*Np*\|*fp* − ∑*q*∈*Np* *fq* = ∑ *q*∈*Np* *vpq*. (8)
>
> Equations (7) form a classical, sparse (banded), symmetric,
> positive-definite system. Because of the arbitrary shape of bound-ary
> ∂Ω, we must use well-known iterative solvers. Results shown in this
> paper have been computed using either Gauss-Seidel iteration with
> successive overrelaxation or V-cycle multigrid. Both methods are fast
> enough for interactive editing of medium size color image regions,
> e.g., 0.4 s. per system on a Pentium 4 for a disk-shaped re-gion of
> 60,000 pixels. As demonstrated in \[Bolz et al. 2003\], multi-grid
> implementation on a GPU will provide a solution for much larger
> regions.
>
> 3 Seamless cloning
>
> Importing gradients The basic choice for the guidance field **v** is a
> gradient field taken directly from a source image. Denoting by *g*
> this source image, the interpolation is performed under the guidance
> of\
> **v** = ∇*g*, (9)
>
> and (4) now reads
>
> ∆*f* = ∆*g* over Ω, with *f*\|∂Ω = *f*∗\|∂Ω. (10)
>
> 315

![](vertopal_952755e54fbc4a5886c9aadf32ad9f47/media/image14.png){width="3.236111111111111in"
height="1.3472222222222223in"}![](vertopal_952755e54fbc4a5886c9aadf32ad9f47/media/image15.png){width="0.8041655730533683in"
height="1.3240682414698162in"}![](vertopal_952755e54fbc4a5886c9aadf32ad9f47/media/image16.png){width="0.7625in"
height="1.3213779527559055in"}![](vertopal_952755e54fbc4a5886c9aadf32ad9f47/media/image17.png){width="0.7625in"
height="1.3213779527559055in"}![](vertopal_952755e54fbc4a5886c9aadf32ad9f47/media/image18.png){width="0.7625in"
height="1.3213779527559055in"}![](vertopal_952755e54fbc4a5886c9aadf32ad9f47/media/image19.png){width="1.2444444444444445in"
height="1.3912817147856518in"}![](vertopal_952755e54fbc4a5886c9aadf32ad9f47/media/image20.png){width="1.245832239720035in"
height="1.392833552055993in"}![](vertopal_952755e54fbc4a5886c9aadf32ad9f47/media/image21.png){width="3.263888888888889in"
height="1.5555555555555556in"}![](vertopal_952755e54fbc4a5886c9aadf32ad9f47/media/image22.png){width="0.6666666666666666in"
height="0.3658814523184602in"}![](vertopal_952755e54fbc4a5886c9aadf32ad9f47/media/image23.png){width="0.6652777777777777in"
height="0.3759514435695538in"}![](vertopal_952755e54fbc4a5886c9aadf32ad9f47/media/image24.png){width="0.6611100174978127in"
height="0.7391174540682415in"}![](vertopal_952755e54fbc4a5886c9aadf32ad9f47/media/image25.png){width="0.6611111111111111in"
height="0.8763560804899387in"}![](vertopal_952755e54fbc4a5886c9aadf32ad9f47/media/image26.png){width="0.6583333333333333in"
height="0.9754363517060367in"}![](vertopal_952755e54fbc4a5886c9aadf32ad9f47/media/image27.png){width="1.2416666666666667in"
height="1.839747375328084in"}![](vertopal_952755e54fbc4a5886c9aadf32ad9f47/media/image28.png){width="1.2416655730533683in"
height="1.8397451881014872in"}![](vertopal_952755e54fbc4a5886c9aadf32ad9f47/media/image29.png){width="3.25in"
height="1.9027777777777777in"}

+--------+--------+--------+--------+--------+--------+--------+--------+
| s      | destin | c      | se     | ![     | ![     | ![     | ![     |
| ources | ations | loning | amless | ](vert | ](vert | ](vert | ](vert |
|        |        |        | c      | opal_9 | opal_9 | opal_9 | opal_9 |
|        |        |        | loning | 52755e | 52755e | 52755e | 52755e |
|        |        |        |        | 54fbc4 | 54fbc4 | 54fbc4 | 54fbc4 |
|        |        |        |        | a5886c | a5886c | a5886c | a5886c |
|        |        |        |        | 9aadf3 | 9aadf3 | 9aadf3 | 9aadf3 |
|        |        |        |        | 2ad9f4 | 2ad9f4 | 2ad9f4 | 2ad9f4 |
|        |        |        |        | 7/medi | 7/medi | 7/medi | 7/medi |
|        |        |        |        | a/imag | a/imag | a/imag | a/imag |
|        |        |        |        | e2.png | e3.png | e4.png | e5.png |
|        |        |        |        | ){widt | ){widt | ){widt | ){widt |
|        |        |        |        | h="0.6 | h="0.6 | h="0.9 | h="0.9 |
|        |        |        |        | 222211 | 180555 | 291666 | 277766 |
|        |        |        |        | 286089 | 555555 | 666666 | 841644 |
|        |        |        |        | 239in" | 556in" | 667in" | 794in" |
|        |        |        |        | height | height | height | height |
|        |        |        |        | ="0.66 | ="0.87 | ="1.31 | ="1.31 |
|        |        |        |        | 666666 | 916666 | 944444 | 944444 |
|        |        |        |        | 666666 | 666666 | 444444 | 444444 |
|        |        |        |        | 66in"} | 67in"} | 44in"} | 44in"} |
+========+========+========+========+========+========+========+========+
|        |        |        |        | source |        |        |        |
|        |        |        |        | /desti |        |        |        |
|        |        |        |        | nation |        |        |        |
+--------+--------+--------+--------+--------+--------+--------+--------+
|        |        |        |        | color  |        |        | mono   |
|        |        |        |        | tr     |        |        | chrome |
|        |        |        |        | ansfer |        |        | tr     |
|        |        |        |        |        |        |        | ansfer |
+--------+--------+--------+--------+--------+--------+--------+--------+
|        |        |        |        | >      |        |        |        |
|        |        |        |        | Figure |        |        |        |
|        |        |        |        | > 5:   |        |        |        |
|        |        |        |        | >      |        |        |        |
|        |        |        |        | **Mono |        |        |        |
|        |        |        |        | chrome |        |        |        |
|        |        |        |        | >      |        |        |        |
|        |        |        |        |  trans |        |        |        |
|        |        |        |        | fer**. |        |        |        |
|        |        |        |        | > In   |        |        |        |
|        |        |        |        | > some |        |        |        |
|        |        |        |        | >      |        |        |        |
|        |        |        |        | cases, |        |        |        |
|        |        |        |        | > such |        |        |        |
|        |        |        |        | > as   |        |        |        |
|        |        |        |        | > tex- |        |        |        |
|        |        |        |        | >      |        |        |        |
|        |        |        |        | > ture |        |        |        |
|        |        |        |        | > tra  |        |        |        |
|        |        |        |        | nsfer, |        |        |        |
|        |        |        |        | > the  |        |        |        |
|        |        |        |        | > part |        |        |        |
|        |        |        |        | > of   |        |        |        |
|        |        |        |        | > the  |        |        |        |
|        |        |        |        | >      |        |        |        |
|        |        |        |        | source |        |        |        |
|        |        |        |        | >      |        |        |        |
|        |        |        |        |  color |        |        |        |
|        |        |        |        | > rem  |        |        |        |
|        |        |        |        | aining |        |        |        |
|        |        |        |        | >      |        |        |        |
|        |        |        |        |  after |        |        |        |
|        |        |        |        | > se   |        |        |        |
|        |        |        |        | amless |        |        |        |
+--------+--------+--------+--------+--------+--------+--------+--------+
|        |        |        |        | > c    |        |        |        |
|        |        |        |        | loning |        |        |        |
|        |        |        |        | >      |        |        |        |
|        |        |        |        |  might |        |        |        |
|        |        |        |        | > be   |        |        |        |
|        |        |        |        | >      |        |        |        |
|        |        |        |        | undesi |        |        |        |
|        |        |        |        | rable. |        |        |        |
|        |        |        |        | > This |        |        |        |
|        |        |        |        | > is   |        |        |        |
|        |        |        |        | >      |        |        |        |
|        |        |        |        |  fixed |        |        |        |
|        |        |        |        | > by   |        |        |        |
|        |        |        |        | > t    |        |        |        |
|        |        |        |        | urning |        |        |        |
|        |        |        |        | > the  |        |        |        |
|        |        |        |        | >      |        |        |        |
|        |        |        |        | source |        |        |        |
+--------+--------+--------+--------+--------+--------+--------+--------+
|        |        |        |        | >      |        |        |        |
|        |        |        |        |  image |        |        |        |
|        |        |        |        | > mono |        |        |        |
|        |        |        |        | chrome |        |        |        |
|        |        |        |        | >      |        |        |        |
|        |        |        |        |  befor |        |        |        |
|        |        |        |        | ehand. |        |        |        |
+--------+--------+--------+--------+--------+--------+--------+--------+

Figs. 6 and 7.

+-----------+-----------+-----------+-----------+-----------+-----------+
| > so      | cloning   | seamless  | ![        | ![        | > ![      |
| urces/des |           | cloning   | ](vertopa | ](vertopa | ](vertopa |
| tinations |           |           | l_952755e | l_952755e | l_952755e |
|           |           |           | 54fbc4a58 | 54fbc4a58 | 54fbc4a58 |
|           |           |           | 86c9aadf3 | 86c9aadf3 | 86c9aadf3 |
|           |           |           | 2ad9f47/m | 2ad9f47/m | 2ad9f47/m |
|           |           |           | edia/imag | edia/imag | edia/imag |
|           |           |           | e6.png){w | e7.png){w | e8.png){w |
|           |           |           | idth="0.5 | idth="1.0 | idth="1.0 |
|           |           |           | 861111111 | 916666666 | 916666666 |
|           |           |           | 111111in" | 666666in" | 666666in" |
|           |           |           | hei       | hei       | > hei     |
|           |           |           | ght="0.54 | ght="1.01 | ght="1.01 |
|           |           |           | 444444444 | 388888888 | 388888888 |
|           |           |           | 44444in"} | 88888in"} | 88888in"} |
+===========+===========+===========+===========+===========+===========+
| Figure 3: |           |           |           |           |           |
| **Ins     |           |           |           |           |           |
| ertion**. |           |           |           |           |           |
| The power |           |           |           |           |           |
| of the    |           |           |           |           |           |
| method is |           |           |           |           |           |
| fully     |           |           |           |           |           |
| expressed |           |           |           |           |           |
| when      |           |           |           |           |           |
| inserting |           |           |           |           |           |
| objects   |           |           |           |           |           |
| with      |           |           |           |           |           |
| complex   |           |           |           |           |           |
| outlines  |           |           |           |           |           |
| into a    |           |           |           |           |           |
| new back- |           |           |           |           |           |
+-----------+-----------+-----------+-----------+-----------+-----------+
| ground.   |           |           |           |           |           |
| Because   |           |           |           |           |           |
| of the    |           |           |           |           |           |
| drastic   |           |           |           |           |           |
| di        |           |           |           |           |           |
| fferences |           |           |           |           |           |
| between   |           |           |           |           |           |
| the       |           |           |           |           |           |
| source    |           |           |           |           |           |
| and       |           |           |           |           |           |
+-----------+-----------+-----------+-----------+-----------+-----------+
| the       |           |           |           |           |           |
| des       |           |           |           |           |           |
| tination, |           |           |           |           |           |
| standard  |           |           |           |           |           |
| image     |           |           |           |           |           |
| cloning   |           |           |           |           |           |
| cannot be |           |           |           |           |           |
| used in   |           |           |           |           |           |
| this      |           |           |           |           |           |
| case.     |           |           |           |           |           |
+-----------+-----------+-----------+-----------+-----------+-----------+
|           |           |           |           | \(a\)     | > \(b\)   |
|           |           |           |           | co        | >         |
|           |           |           |           | lor-based |  seamless |
|           |           |           |           | cutout    | > cloning |
|           |           |           |           | and paste |           |
+-----------+-----------+-----------+-----------+-----------+-----------+
|           |           |           | ![        | ![]       | > ![]     |
|           |           |           | ](vertopa | (vertopal | (vertopal |
|           |           |           | l_952755e | _952755e5 | _952755e5 |
|           |           |           | 54fbc4a58 | 4fbc4a588 | 4fbc4a588 |
|           |           |           | 86c9aadf3 | 6c9aadf32 | 6c9aadf32 |
|           |           |           | 2ad9f47/m | ad9f47/me | ad9f47/me |
|           |           |           | edia/imag | dia/image | dia/image |
|           |           |           | e9.png){w | 10.png){w | 11.png){w |
|           |           |           | idth="0.5 | idth="1.0 | idth="1.0 |
|           |           |           | 847222222 | 916666666 | 916666666 |
|           |           |           | 222223in" | 666666in" | 666666in" |
|           |           |           | hei       | height="1 | >         |
|           |           |           | ght="0.54 | .0125in"} | height="1 |
|           |           |           | 305446194 |           | .0125in"} |
|           |           |           | 22573in"} |           |           |
+-----------+-----------+-----------+-----------+-----------+-----------+
|           |           |           |           | \(c\)     | > \(d\)   |
|           |           |           |           | seamless  | > mixed   |
|           |           |           |           | cloning   | >         |
|           |           |           |           | and       |  seamless |
|           |           |           |           | de        | > cloning |
|           |           |           |           | stination |           |
|           |           |           |           | av-       |           |
+-----------+-----------+-----------+-----------+-----------+-----------+
|           |           |           |           | > eraged  |           |
+-----------+-----------+-----------+-----------+-----------+-----------+

> Figure 6: **Inserting objects with holes**. (a) The classic method,\
> color-based selection and alpha masking might be time consuming\
> and often leaves an undesirable halo; (b-c) seamless cloning, even\
> averaged with the original image, is not effective; (d) mixed seam-\
> less cloning based on a loose selection proves effective.

+-----------------------+-----------------------+-----------------------+
| > source/destination  | cloning               | > seamless cloning    |
+=======================+=======================+=======================+
+-----------------------+-----------------------+-----------------------+

![](vertopal_952755e54fbc4a5886c9aadf32ad9f47/media/image12.png){width="0.7791655730533683in"
height="0.5847222222222223in"}![](vertopal_952755e54fbc4a5886c9aadf32ad9f47/media/image13.png){width="1.95in"
height="1.4625in"}

swapped textures

Figure 4: **Feature exchange**. Seamless cloning allows the user to
replace easily certain features of one object by alternative features.
In the second example of texture swapping multiple broad strokes (not
shown) were used.

+-----------------+-----------------+-----------------+-----------------+
| The discrete    |                 |                 | \(13\)          |
| counterpart of  |                 |                 |                 |
| this guidance   |                 |                 |                 |
| field is:       |                 |                 |                 |
+=================+=================+=================+=================+
| *vpq* =         | >   *f* ∗*p* −  | > if \|*f*∗*p*  |                 |
|                 | > *f* ∗*q*      | > − *f*∗*q* \|  |                 |
|                 |                 | > \> \|*gp*     |                 |
|                 |                 | > −*gq*\|,      |                 |
|                 |                 | > otherwise,    |                 |
+-----------------+-----------------+-----------------+-----------------+

for all ⟨*p*,*q*⟩. The effect of this guidance field is demonstrated in

> 316

![](vertopal_952755e54fbc4a5886c9aadf32ad9f47/media/image32.png){width="3.2916666666666665in"
height="1.375in"}![](vertopal_952755e54fbc4a5886c9aadf32ad9f47/media/image33.png){width="0.7055555555555556in"
height="0.5885454943132109in"}![](vertopal_952755e54fbc4a5886c9aadf32ad9f47/media/image34.png){width="0.7013888888888888in"
height="0.75040135608049in"}![](vertopal_952755e54fbc4a5886c9aadf32ad9f47/media/image35.png){width="1.2472222222222222in"
height="1.3343777340332459in"}![](vertopal_952755e54fbc4a5886c9aadf32ad9f47/media/image36.png){width="1.2472222222222222in"
height="1.3343777340332459in"}![](vertopal_952755e54fbc4a5886c9aadf32ad9f47/media/image37.png){width="1.1666666666666667in"
height="1.1584317585301838in"}![](vertopal_952755e54fbc4a5886c9aadf32ad9f47/media/image38.png){width="2.4027777777777777in"
height="2.5416666666666665in"}![](vertopal_952755e54fbc4a5886c9aadf32ad9f47/media/image39.png){width="1.1638877952755906in"
height="1.3415201224846893in"}![](vertopal_952755e54fbc4a5886c9aadf32ad9f47/media/image40.png){width="1.163888888888889in"
height="1.3415212160979877in"}![](vertopal_952755e54fbc4a5886c9aadf32ad9f47/media/image41.png){width="1.1666666666666667in"
height="1.1584317585301838in"}

+-----------------------+-----------------------+-----------------------+
| > source/destination  | seamless cloning      | mixed seamless        |
|                       |                       | cloning               |
+=======================+=======================+=======================+
+-----------------------+-----------------------+-----------------------+

Figure 8: **Inserting one object close to another**. With seamless
cloning, an object in the destination image touching the selected region
Ω bleeds into it. Bleeding is inhibited by using mixed gradi-ents as the
guidance field.

4 Selection editing

In the two previous sections, the guidance field depended, partly or
wholly, on the gradient field of a source image *g*. Alternatively,
in-place image transformations can be defined by using a guidance field
depending entirely on the original image. Based on this idea, this
section details texture flattening, spatially selective illumination
changes, background or foreground color modifications, and seam-less
tiling. The first two effects rely on non-linear modifications of the
original gradient field ∇*f*∗in the selected region. The latter effects
rely on in-place seamless cloning after the original image has been
modified either inside the domain, providing a new source image, or
outside, providing new boundary conditions.

+-----------------------------------+-----------------------------------+
| ![](vertopal_952755e54fbc4a       | > ![](vertopal_952755e54fbc4a     |
| 5886c9aadf32ad9f47/media/image30. | 5886c9aadf32ad9f47/media/image31. |
| png){width="1.2805555555555554in" | png){width="1.2805544619422573in" |
| height="1.6694444444444445in"}    | > height="1.6694444444444445in"}  |
+===================================+===================================+
| > Figure 9: **Texture             |                                   |
| > flattening**. By retaining only |                                   |
| > the gradients at edge           |                                   |
| > locations, before integrating   |                                   |
| > with the Poisson solver, one    |                                   |
| > washes out the texture of the   |                                   |
| > selected region, giving its     |                                   |
| > contents a flat aspect.         |                                   |
+-----------------------------------+-----------------------------------+

> A natural extension is to restrict the correction to a selected
> re-gion Ω, using appropriate Dirichlet conditions on ∂Ω. Using a
> sim-plified version of the Fattal *et al.* transformation \[Fattal et
> al. 2002\], the guidance field is defined in the log-domain by:
>
> **v** = αβ\|∇*f*∗\|−β∇*f*∗, (16)
>
> with α = 0.2 times the average gradient norm of *f*∗over Ω, and β =
> 0.2. As shown in Fig. 10, this tool can be used for instance to
> correct an under-exposed object of interest, or to reduce specular
> reflections.

Texture flattening The image gradient ∇*f*∗is passed through a

sparse sieve that retains only the most salient features:

> for all **x** ∈ Ω, **v**(**x**) = *M*(**x**)∇*f*∗(**x**), (14)

where *M* is a binary mask turned on at a few locations of interest.

> A good choice for *M* is an edge detector, in which case the dis-

crete version of (14), to be plugged into (7), is:

+-------------+-------------+-------------+-------------+-------------+
| *vpq* =     |             | 0 *fp* −    | if an edge  | > \(15\)    |
|             |             | *fq*        | lies        |             |
|             |             |             | between *p* |             |
|             |             |             | and *q*,    |             |
+=============+=============+=============+=============+=============+
|             |             |             | >           |             |
|             |             |             |  otherwise, |             |
+-------------+-------------+-------------+-------------+-------------+

for all ⟨*p*,*q*⟩. As shown in Fig. 9, the content of the selection Ω
gets a flattened appearance, with small grain details washed out, and
the main structure preserved. The extent of this effect depends
obviously on the sparsity of the sieve. The more selective the edge
detector, the sparser the edge map, and the more pronounced the effect.

Note that this instance of Poisson editing has strong connections with
the contour-domain editing system of Elder and Goldberg \[El-der and
Goldberg 2001\]. The difference is that we specify approx-imately the
gradient vectors at edge locations through sparse guid-ance (14),
whereas their system relies on an exact specification of color values on
both sides of each edgel.

Local illumination changes As pointed out by the authors, the method of
\[Fattal et al. 2002\] is not limited to HDR images and can be applied
to ordinary images in order to modify smoothly their dy-namic range.
First, the gradient field of the logarithm of the image is transformed
in order to reduce the large gradients and to increase the small ones.
The transformed vector field **v** is then used to recon-struct the
logarithm of the image, *f*, by solving the Poisson equa-tion ∆*f* =
div**v** over the whole image domain under the Neumann boundary
conditions.

> 317

![](vertopal_952755e54fbc4a5886c9aadf32ad9f47/media/image44.png){width="1.0472222222222223in"
height="1.5572692475940508in"}![](vertopal_952755e54fbc4a5886c9aadf32ad9f47/media/image45.png){width="1.0472222222222223in"
height="1.5572692475940508in"}![](vertopal_952755e54fbc4a5886c9aadf32ad9f47/media/image46.png){width="1.0472222222222223in"
height="1.5572692475940508in"}![](vertopal_952755e54fbc4a5886c9aadf32ad9f47/media/image47.png){width="3.25in"
height="1.5694444444444444in"}

Figure 11: **Local color changes**. Left: original image showing
selection Ω surrounding loosely an object of interest; center:
back-ground decolorization done by setting *g* to the original color
image and *f*∗to the luminance of *g*; right: recoloring the object of
interest by multiplying the RGB channels of the original image by 1.5,
0.5, and 0.5 respectively to form the source image.

> ![](vertopal_952755e54fbc4a5886c9aadf32ad9f47/media/image42.png){width="1.6041655730533684in"
> height="1.0694444444444444in"}![](vertopal_952755e54fbc4a5886c9aadf32ad9f47/media/image43.png){width="1.6055555555555556in"
> height="1.0694444444444444in"}

Figure 12: **Seamless tiling**. Setting periodic boundary values on the
border of a rectangular region before integrating with the Poisson
solver yields a tileable image.

performed by precisely selecting an object and then setting its
com-plement to monochrome. In contrast, Poisson editing frees the user
from the tedium of precise selection: given a source color image *g*,
(a) the destination image *f*∗is set to be the luminance channel from
*g*, (b) the user selects a region Ω containing the object, and this may
be somewhat bigger than the actual object, and (c) the Pois-son equation
(10) is solved in each color channel. An example is presented in Fig.
11. Note that, although the result seems to offer also a precise
segmentation of the object for free, this is not actually the case as
there is some residual contamination of the destination image outside
the object.

Conversely Poisson image editing can be used to modify the color of a
loosely selected object. Before solving the Poisson equa-tion (10), the
original image is copied to the destination *f*∗and a version with
modified colors is copied to the source *g*, see Fig. 11.

Seamless tiling When the domain Ω is rectangular, its content can be
made tileable by enforcing periodic boundary conditions with the Poisson
solver. The source image *g* is the original im-age, and the boundary
conditions are derived from the boundary values of *g*, such that
opposite sides of the rectangular domain cor-respond to identical
Dirichlet conditions. In Fig. 12, we have cho-sen *f*∗\
north= *f* ∗south= 0.5(*g*north + *g*south), and similarly for the east
and west borders.

> image inside the selection, such as texture, illumination, or color.
> An important common characteristic of all these tools is that there is
> no need for precise object delineation, in contrast with the classic
> tools that address similar tasks. This is a valuable feature, whether
> one is interested in small touch-up operations or in complex
> photo-montages.
>
> Although not illustrated in this paper, it is clear that the cloning
> facilities described in Section 3 can be combined with the editing
> ones introduced in Section 4. It is for instance possible to insert an
> object while flattening its texture to make it match the style of a
> texture-free destination.
>
> Finally, it is worth noting that the range of editing facilities
> de-rived in this paper from the same generic framework could prob-ably
> be extended further. Appearance changes could for instance also deal
> with the sharpness of objects of interest, thus allowing the user to
> make apparent changes of focus.
>
> Image credits Two landscapes and swimming bear in Fig. 3, flower in
> Fig.11: from Corel Professional Photos, copyright c⃝2003 Microsoft
> Research and its licensors, all rights reserved; rainbow in Fig. 7
> courtesy Professor James B. Kaler, University of Illinois.
>
> References
>
> ADOBEc⃝. 2002. *Photoshop* c⃝ *7.0 User Guide*. Adobe Systems Incorpo-
> rated.
>
> BALLESTER, C., BERTALMIO, M., CASELLES, V., SAPIRO, G., AND VERDERA,
> J. 2001. Filling-in by Joint Interpolation of Vector Fields and Gray
> Levels. *IEEE Trans. Image Processing 10*, 8, 1200--1211.
>
> BARRET, A., AND CHENEY, A. 2002. Object-Based Image Editing. *ACM*
> *Transactions on Graphics 21*, 3, 777--784.
>
> BERTALMIO, M., SAPIRO, G., CASELLES, V., AND BALLESTER, C. 2000. Image
> Inpainting. In *Proceedings of ACM SIGGRAPH 2000*, ACM Press / ACM
> SIGGRAPH, New-York, E. Fiume, Ed., Computer Graphics Proceedings,
> Annual Conference Series, ACM, 417--424.
>
> BOLZ, J., FARMER, I., GRINPSUN, E., AND SCHR ¨ODER, P. 2003. Sparse
> Matrix Solvers on the GPU: Conjugate Gradients and Multigrid. *ACM
> Transactions on Graphics*. to appear.
>
> BORNARD, R., LECAN, E., LABORELLI, L., AND CHENOT, J.-H. 2002. Missing
> Data Correction in Still Images and Image Sequences. In *Proc. ACM
> International Conference on Multimedia*.
>
> BURT, P., AND ADELSON, E. 1983. A Multiresolution Spline with
> Appli-cation to Image Mosaics. *ACM Transactions on Graphics 2*, 4,
> 217--236.
>
> EFROS, A., AND FREEMAN, W. 2001. Image Quilting for Texture Synthe-sis
> and Transfer. In *Proceedings of ACM SIGGRAPH 2001*, ACM Press / ACM
> SIGGRAPH, New-York, E. Fiume, Ed., Computer Graphics Pro-ceedings,
> Annual Conference Series, ACM, 341--346.
>
> EFROS, A., AND LEUNG, T. 1999. Texture Synthesis by Non-Parametric
> Sampling. In *Proc. Int. Conf. Computer Vision*, 1033--1038.
>
> ELDER, J., AND GOLDBERG, R. 2001. Image Editing in the Contour Domain.
> *IEEE Trans. Pattern Anal. Machine Intell. 23*, 3, 291--296.
>
> FATTAL, R., LISCHINSKI, D., AND WERMAN, M. 2002. Gradient Domain High
> Dynamic Range Compression. *ACM Transactions on Graphics 21*, 3,
> 249--256.
>
> HERTZMANN, A., JACOBS, C., OLIVER, N., CURLESS, B., AND SALESIN, D.
> 2001. Image Analogies. In *Proceedings of SIGGRAPH* *2001*, ACM Press
> / ACM SIGGRAPH, New-York, E. Fiume, Ed., Com- puter Graphics
> Proceedings, Annual Conference Series, ACM, 327--340.

+-----------------------+-----------------------+-----------------------+
| 5                     | > Conclusion          | LAND, E., AND MCCANN, |
|                       |                       | J. 1971. Ligthness    |
|                       |                       | and Retinex Theory.   |
|                       |                       | *J. Opt.*             |
+=======================+=======================+=======================+
|                       |                       | *Soc. Amer. 61*,      |
|                       |                       | 1--11.                |
+-----------------------+-----------------------+-----------------------+

Using the generic framework of guided interpolation, we have in-troduced
a variety of tools to edit in a seamless and effortless manner the
contents of an image selection. The extent of possi-ble changes ranges
from replacement by, or mixing with, another source image region, to
alterations of some aspects of the original

> 318
