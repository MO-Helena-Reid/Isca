! This file was translated by Helena Reid
! as part of the initial release of the 
! PAPILLON stochastic physics scheme,
! under the BSD 3-clause license,
! which should have been provided with this file.
!
! A description of the scheme may be found in preprint here:
! https://doi.org/10.5194/egusphere-2025-6312
!
! Implementation of simplex noise in 4d
! Translated to python then to Fortran 90, adapted from DOI:10.13140/RG.2.1.3369.6488
 
module simplex_noise_mod
    implicit none
    contains
   function snoise4d(v) result(noise)
        implicit none
        real :: noise
        real, dimension(4), intent(in) :: v
        integer, parameter, dimension(64, 4) :: SIMPLEX = reshape([ &
        0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, &
        0, 0, 1, 0, 1, 0, 0, 0, 2, 2, 1, 1, 0, 0, 0, 2, 0, 2, 0, 0, 0, 0, &
        0, 0, 0, 0, 2, 0, 0, 0, 3, 3, 0, 3, 2, 0, 0, 0, 3, 0, 3, 3, &
        1, 1, 0, 2, 0, 0, 0, 2, 2, 0, 3, 3, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, &
        0, 0, 2, 0, 3, 0, 0, 0, 3, 3, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 2, 2, &
        2, 3, 0, 3, 0, 0, 0, 3, 1, 0, 1, 2, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, &
        0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 3, 0, 0, 0, 3, 0, 3, 0, 0, 0, 0, &
        0, 0, 0, 0, 1, 0, 0, 0, 1, 2, 0, 2, 0, 0, 0, 0, 0, 0, 0, 1, &
        3, 2, 0, 1, 0, 0, 0, 0, 3, 0, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
        0, 0, 3, 0, 2, 0, 0, 0, 1, 0, 3, 2, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, &
        0, 0, 0, 0, 3, 0, 0, 0, 2, 1, 0, 0, 3, 0, 0, 0, 2, 0, 1, 0 &
        ], [64, 4])
        integer, parameter, dimension(32, 4) :: GRAD4 = reshape([ &
            0,  0,  0,  0,  0,  0,  0,  0,  1,  1,  1,  1, -1, -1, -1, -1,&
            1,  1,  1,  1, -1, -1, -1, -1,  1,  1,  1,  1, -1, -1, -1, -1,&
            1,  1,  1,  1, -1, -1, -1, -1,  0,  0,  0,  0,  0,  0,  0,  0,&
            1,  1, -1, -1,  1,  1, -1, -1,  1,  1, -1, -1,  1,  1, -1, -1,&
            1,  1, -1, -1,  1,  1, -1, -1,  1,  1, -1, -1,  1,  1, -1, -1,&
            0,  0,  0,  0,  0,  0,  0,  0,  1, -1,  1, -1,  1, -1,  1, -1,&
            1, -1,  1, -1,  1, -1,  1, -1,  1, -1,  1, -1,  1, -1,  1, -1,&
            1, -1,  1, -1,  1, -1,  1, -1,  0,  0,  0,  0,  0,  0,  0,  0 &
        ], [32, 4])
        integer, parameter :: p(256) = [ &
            151, 160, 137, 91, 90, 15, 131, 13, 201, 95, 96, 53, 194, 233, 7, 225, &
            140, 36, 103, 30, 69, 142, 8, 99, 37, 240, 21, 10, 23, 190, 6, 148, &
            247, 120, 234, 75, 0, 26, 197, 62, 94, 252, 219, 203, 117, 35, 11, 32, &
            57, 177, 33, 88, 237, 149, 56, 87, 174, 20, 125, 136, 171, 168, 68, 175, &
            74, 165, 71, 134, 139, 48, 27, 166, 77, 146, 158, 231, 83, 111, 229, 122, &
            60, 211, 133, 230, 220, 105, 92, 41, 55, 46, 245, 40, 244, 102, 143, 54, &
            65, 25, 63, 161, 1, 216, 80, 73, 209, 76, 132, 187, 208, 89, 18, 169, &
            200, 196, 135, 130, 116, 188, 159, 86, 164, 100, 109, 198, 173, 186, 3, 64, &
            52, 217, 226, 250, 124, 123, 5, 202, 38, 147, 118, 126, 255, 82, 85, 212, &
            207, 206, 59, 227, 47, 16, 58, 17, 182, 189, 28, 42, 223, 183, 170, 213, &
            119, 248, 152, 2, 44, 154, 163, 70, 221, 153, 101, 155, 167, 43, 172, 9, &
            129, 22, 39, 253, 19, 98, 108, 110, 79, 113, 224, 232, 178, 185, 112, 104, &
            218, 246, 97, 228, 251, 34, 242, 193, 238, 210, 144, 12, 191, 179, 162, 241, &
            81, 51, 145, 235, 249, 14, 239, 107, 49, 192, 214, 31, 181, 199, 106, 157, &
            184, 84, 204, 176, 115, 121, 50, 45, 127, 4, 150, 254, 138, 236, 205, 93, &
            222, 114, 67, 29, 24, 72, 243, 141, 128, 195, 78, 66, 215, 61, 156, 180 &
        ]
        real, parameter :: F4 = (sqrt(5.0) - 1.0) / 4.0
        real, parameter :: G4 = (5.0 - sqrt(5.0)) / 20.0
        integer :: i, c1, c2, c3, c4, c5, c6, c, gi0, gi1, gi2, gi3, gi4
        integer, dimension(4) :: hc_origin, corner1, corner2, corner3, &
                                 origin_mod_256
        integer, parameter :: PERM(512) = [ (p(i), i=1, 256), (p(i), i=1, 256) ]
        real :: n0, n1, n2, n3, n4, &
                           s, t, t0, t1, t2, t3, t4
        real, dimension(4) :: v0, v1, v2, v3, v4
        s = sum(v) * F4
        hc_origin(:) = floor(v(:) + s)

        t = sum(hc_origin) * G4
        
        v0 = v - hc_origin + t
        
        ! We're trying to find which simplex in the hypercube we're in.
        ! As there are 64 possible orderings of four numbers
        ! in a 4-vector, we need a 6-bit integer to traverse them.
        c1 = merge(32, 0, v0(1) > v0(2))
        c2 = merge(16, 0, v0(1) > v0(3))
        c3 = merge(8, 0, v0(2) > v0(3))
        c4 = merge(4, 0, v0(1) > v0(4))
        c5 = merge(2, 0, v0(2) > v0(4))
        c6 = merge(1, 0, v0(3) > v0(4))
        c = c1 + c2 + c3 + c4 + c5 + c6 + 1
        corner1(:) = merge(1, 0, SIMPLEX(c, :) >= 3)
        corner2(:) = merge(1, 0, SIMPLEX(c, :) >= 2)
        corner3(:) = merge(1, 0, SIMPLEX(c, :) >= 1)
        v1 = v0 - corner1 + G4
        v2 = v0 - corner2 + 2.0 * G4
        v3 = v0 - corner3 + 3.0 * G4
        v4 = v0 - 1.0 + 4.0 * G4
        
        origin_mod_256 = modulo(hc_origin, 256) + 1
        
        ! Now we compute the hashed gradient indices of the five simplex corners
        gi0 = modulo(PERM(origin_mod_256(1) + PERM(origin_mod_256(2) + &
                   PERM(origin_mod_256(3) + PERM(origin_mod_256(4))))), 32)
        ! Getting crashes in PERM() below due to SIGSEGV invalid memory reference.
        ! Let's check if we're passing out of bound indices to PERM
!        if (origin_mod_256(4) + corner1(4) > 512) then
!            print*, 'origin_mod_256(4) + corner1(4) > 512'
!            print*, 'origin_mod_256(4) + corner1(4)', origin_mod_256(4) + corner1(4)
!            print*, 'origin_mod_256(4)', origin_mod_256(4)
!            print*, 'corner1(4)', corner1(4)
!            print*, 'hc_origin(4)', hc_origin(4)
!        end if
!        if (origin_mod_256(4) + corner1(4) < 1) then
!            print*, 'origin_mod_256(4) + corner1(4) < 1'
!            print*, 'origin_mod_256(4) + corner1(4)', origin_mod_256(4) + corner1(4)
!            print*, 'origin_mod_256(4)', origin_mod_256(4)
!            print*, 'corner1(4)', corner1(4)
!            print*, 'hc_origin(4)', hc_origin(4)
!        end if
!        if (origin_mod_256(4) + corner2(4) > 512) then
!            print*, 'origin_mod_256(4) + corner2(4) > 512'
!            print*, 'origin_mod_256(4) + corner2(4)', origin_mod_256(3) + corner2(3)
!            print*, 'origin_mod_256(4)', origin_mod_256(4)
!            print*, 'corner2(4)', corner2(4)
!            print*, 'hc_origin(4)', hc_origin(4)
!        end if
!        if (origin_mod_256(4) + corner2(4) < 1) then
!            print*, 'origin_mod_256(4) + corner2(4) < 1'
!            print*, 'origin_mod_256(4) + corner2(4)', origin_mod_256(4) + corner2(4)
!            print*, 'origin_mod_256(4)', origin_mod_256(4)
!            print*, 'corner2(4)', corner2(4)
!            print*, 'hc_origin(4)', hc_origin(4)
!        end if
!        if (origin_mod_256(4) + corner3(4) > 512) then
!            print*, 'origin_mod_256(4) + corner3(4) > 512'
!            print*, 'origin_mod_256(4) + corner3(4)', origin_mod_256(4) + corner3(4)
!            print*, 'origin_mod_256(4)', origin_mod_256(4)
!            print*, 'corner3(4)', corner3(4)
!            print*, 'hc_origin(4)', hc_origin(4)
!        end if
!        if (origin_mod_256(4) + corner3(4) < 1) then
!            print*, 'origin_mod_256(4) + corner3(4) < 1'
!            print*, 'origin_mod_256(4) + corner3(4)', origin_mod_256(4) + corner3(4)
!            print*, 'origin_mod_256(4)', origin_mod_256(4)
!            print*, 'corner3(4)', corner3(4)
!        end if
!        if (origin_mod_256(3) + corner1(3) + PERM(origin_mod_256(4) + corner1(4)) > 512) then
!            print*, 'origin_mod_256(3) + corner1(3) + PERM(origin_mod_256(4) + corner1(4)) > 512'
!            print*, 'origin_mod_256(3) + corner1(3) + PERM(origin_mod_256(4) + corner1(4))',&
!                    origin_mod_256(3) + corner1(3) + PERM(origin_mod_256(4) + corner1(4))
!            print*, 'origin_mod_256(3)', origin_mod_256(3)
!            print*, 'corner1(3)', corner1(3)
!            print*, 'PERM(origin_mod_256(4) + corner1(4))', PERM(origin_mod_256(4) + corner1(4))
!            print*, 'origin_mod_256(4)', origin_mod_256(4)
!            print*, 'corner1(4)', corner1(4)
!        end if
!        if (origin_mod_256(3) + corner1(3) + PERM(origin_mod_256(4) + corner1(4)) < 1) then
!            print*, 'origin_mod_256(3) + corner1(3) + PERM(origin_mod_256(4) + corner1(4)) < 1'
!            print*, 'origin_mod_256(3) + corner1(3) + PERM(origin_mod_256(4) + corner1(4))',&
!                    origin_mod_256(3) + corner1(3) + PERM(origin_mod_256(4) + corner1(4))
!            print*, 'origin_mod_256(3)', origin_mod_256(3)
!            print*, 'corner1(3)', corner1(3)
!            print*, 'PERM(origin_mod_256(4) + corner1(4))', PERM(origin_mod_256(4) + corner1(4))
!            print*, 'origin_mod_256(4)', origin_mod_256(4)
!            print*, 'corner1(4)', corner1(4)
!            print*, 'hc_origin(4)', hc_origin(4)
!            print*, 'hc_origin', hc_origin
!            print*, 'v', v
!            print*, 's', s
!        end if
!        if (origin_mod_256(3) + corner2(3) + PERM(origin_mod_256(4) + corner2(4)) > 512) then
!            print*, 'origin_mod_256(3) + corner2(3) + PERM(origin_mod_256(4) + corner2(4)) > 512'
!            print*, 'origin_mod_256(3) + corner2(3) + PERM(origin_mod_256(4) + corner2(4))',&
!                    origin_mod_256(3) + corner2(3) + PERM(origin_mod_256(4) + corner2(4))
!            print*, 'origin_mod_256(3)', origin_mod_256(3)
!            print*, 'corner2(3)', corner2(3)
!            print*, 'PERM(origin_mod_256(4) + corner2(4))', PERM(origin_mod_256(4) + corner2(4))
!            print*, 'origin_mod_256(4)', origin_mod_256(4)
!            print*, 'corner2(4)', corner2(4)
!        end if
!        if (origin_mod_256(3) + corner2(3) + PERM(origin_mod_256(4) + corner2(4)) < 1) then
!            print*, 'origin_mod_256(3) + corner2(3) + PERM(origin_mod_256(4) + corner2(4)) < 1'
!            print*, 'origin_mod_256(3) + corner2(3) + PERM(origin_mod_256(4) + corner2(4))',&
!                    origin_mod_256(3) + corner2(3) + PERM(origin_mod_256(4) + corner2(4))
!            print*, 'origin_mod_256(3)', origin_mod_256(3)
!            print*, 'corner2(3)', corner2(3)
!            print*, 'PERM(origin_mod_256(4) + corner2(4))', PERM(origin_mod_256(4) + corner2(4))
!            print*, 'origin_mod_256(4)', origin_mod_256(4)
!            print*, 'corner2(4)', corner2(4)
!        end if
!        if (origin_mod_256(3) + corner3(3) + PERM(origin_mod_256(4) + corner3(4)) > 512) then
!            print*, 'origin_mod_256(3) + corner3(3) + PERM(origin_mod_256(4) + corner3(4)) > 512'
!            print*, 'origin_mod_256(3) + corner3(3) + PERM(origin_mod_256(4) + corner3(4))',&
!                    origin_mod_256(3) + corner3(3) + PERM(origin_mod_256(4) + corner3(4))
!            print*, 'origin_mod_256(3)', origin_mod_256(3)
!            print*, 'corner3(3)', corner3(3)
!            print*, 'PERM(origin_mod_256(4) + corner3(4))', PERM(origin_mod_256(4) + corner3(4))
!            print*, 'origin_mod_256(4)', origin_mod_256(4)
!            print*, 'corner3(4)', corner3(4)
!        end if
!        if (origin_mod_256(3) + corner3(3) + PERM(origin_mod_256(4) + corner3(4)) < 1) then
!            print*, 'origin_mod_256(3) + corner3(3) + PERM(origin_mod_256(4) + corner3(4)) < 1'
!            print*, 'origin_mod_256(3) + corner3(3) + PERM(origin_mod_256(4) + corner3(4))',&
!                    origin_mod_256(3) + corner3(3) + PERM(origin_mod_256(4) + corner3(4))
!            print*, 'origin_mod_256(3)', origin_mod_256(3)
!            print*, 'corner3(3)', corner3(3)
!            print*, 'PERM(origin_mod_256(4) + corner3(4))', PERM(origin_mod_256(4) + corner3(4))
!            print*, 'origin_mod_256(4)', origin_mod_256(4)
!            print*, 'corner3(4)', corner3(4)
!        end if
!        if (origin_mod_256(2) + corner1(2) + PERM(origin_mod_256(3) + corner1(3) + &
!            PERM(origin_mod_256(4) + corner1(4))) > 512) then
!            print*, 'origin_mod_256(2) + corner1(2) + PERM(origin_mod_256(3) + corner1(3) +',&
!                    'PERM(origin_mod_256(4) + corner1(4)) > 512'
!            print*, 'origin_mod_256(2) + corner1(2) + PERM(origin_mod_256(3) + corner1(3) +',&
!                    'PERM(origin_mod_256(4) + corner1(4))',&
!                    origin_mod_256(2) + corner1(2) + PERM(origin_mod_256(3) + corner1(3) + &
!                    PERM(origin_mod_256(4) + corner1(4)))
!            print*, 'origin_mod_256(2)', origin_mod_256(2)
!            print*, 'corner1(2)', corner1(2)
!            print*, 'PERM(origin_mod_256(3) + corner1(3) + PERM(origin_mod_256(4) + corner1(4))',&
!                    PERM(origin_mod_256(3) + corner1(3) + PERM(origin_mod_256(4) + corner1(4)))
!            print*, 'origin_mod_256(3)', origin_mod_256(3)
!            print*, 'corner1(3)', corner1(3)
!            print*, 'PERM(origin_mod_256(4) + corner1(4))', PERM(origin_mod_256(4) + corner1(4))
!            print*, 'origin_mod_256(4)', origin_mod_256(4)
!            print*, 'corner1(4)', corner1(4)
!        end if
!        if (origin_mod_256(2) + corner1(2) + PERM(origin_mod_256(3) + corner1(3) + &
!            PERM(origin_mod_256(4) + corner1(4))) < 1) then
!            print*, 'origin_mod_256(2) + corner1(2) + PERM(origin_mod_256(3) + corner1(3) +', &
!                    'PERM(origin_mod_256(4) + corner1(4)) < 1'
!            print*, 'origin_mod_256(2) + corner1(2) + PERM(origin_mod_256(3) + corner1(3) +', &
!                    'PERM(origin_mod_256(4) + corner1(4))',&
!                    origin_mod_256(2) + corner1(2) + PERM(origin_mod_256(3) + corner1(3) + &
!                    PERM(origin_mod_256(4) + corner1(4)))
!            print*, 'origin_mod_256(2)', origin_mod_256(2)
!            print*, 'corner1(2)', corner1(2)
!            print*, 'PERM(origin_mod_256(3) + corner1(3) + PERM(origin_mod_256(4) + corner1(4))',&
!                    PERM(origin_mod_256(3) + corner1(3) + PERM(origin_mod_256(4) + corner1(4)))
!            print*, 'origin_mod_256(3)', origin_mod_256(3)
!            print*, 'corner1(3)', corner1(3)
!            print*, 'PERM(origin_mod_256(4) + corner1(4))', PERM(origin_mod_256(4) + corner1(4))
!            print*, 'origin_mod_256(4)', origin_mod_256(4)
!            print*, 'corner1(4)', corner1(4)
!        end if
!        if (origin_mod_256(2) + corner2(2) + PERM(origin_mod_256(3) + corner2(3) + &
!            PERM(origin_mod_256(4) + corner2(4))) > 512) then
!            print*, 'origin_mod_256(2) + corner2(2) + PERM(origin_mod_256(3) + corner2(3) +', &
!                    'PERM(origin_mod_256(4) + corner2(4)) > 512'
!            print*, 'origin_mod_256(2) + corner2(2) + PERM(origin_mod_256(3) + corner2(3) +', &
!                    'PERM(origin_mod_256(4) + corner2(4))',&
!                    origin_mod_256(2) + corner2(2) + PERM(origin_mod_256(3) + corner2(3) + &
!                    PERM(origin_mod_256(4) + corner2(4)))
!            print*, 'origin_mod_256(2)', origin_mod_256(2)
!            print*, 'corner2(2)', corner2(2)
!            print*, 'PERM(origin_mod_256(3) + corner2(3) + PERM(origin_mod_256(4) + corner2(4))',&
!                    PERM(origin_mod_256(3) + corner2(3) + PERM(origin_mod_256(4) + corner2(4)))
!            print*, 'origin_mod_256(3)', origin_mod_256(3)
!            print*, 'corner2(3)', corner2(3)
!            print*, 'PERM(origin_mod_256(4) + corner2(4))', PERM(origin_mod_256(4) + corner2(4))
!            print*, 'origin_mod_256(4)', origin_mod_256(4)
!            print*, 'corner2(4)', corner2(4)
!        end if
!        if (origin_mod_256(2) + corner2(2) + PERM(origin_mod_256(3) + corner2(3) + &
!            PERM(origin_mod_256(4) + corner2(4))) < 1) then
!            print*, 'origin_mod_256(2) + corner2(2) + PERM(origin_mod_256(3) + corner2(3) +', &
!                    'PERM(origin_mod_256(4) + corner2(4)) < 1'
!            print*, 'origin_mod_256(2) + corner2(2) + PERM(origin_mod_256(3) + corner2(3) +', &
!                    'PERM(origin_mod_256(4) + corner2(4))',&
!                    origin_mod_256(2) + corner2(2) + PERM(origin_mod_256(3) + corner2(3) + &
!                    PERM(origin_mod_256(4) + corner2(4)))
!            print*, 'origin_mod_256(2)', origin_mod_256(2)
!            print*, 'corner2(2)', corner2(2)
!            print*, 'PERM(origin_mod_256(3) + corner2(3) + PERM(origin_mod_256(4) + corner2(4))',&
!                    PERM(origin_mod_256(3) + corner2(3) + PERM(origin_mod_256(4) + corner2(4)))
!            print*, 'origin_mod_256(3)', origin_mod_256(3)
!            print*, 'corner2(3)', corner2(3)
!            print*, 'PERM(origin_mod_256(4) + corner2(4))', PERM(origin_mod_256(4) + corner2(4))
!            print*, 'origin_mod_256(4)', origin_mod_256(4)
!            print*, 'corner2(4)', corner2(4)
!        end if
!        if (origin_mod_256(2) + corner3(2) + PERM(origin_mod_256(3) + corner3(3) + &
!            PERM(origin_mod_256(4) + corner3(4))) > 512) then
!            print*, 'origin_mod_256(2) + corner3(2) + PERM(origin_mod_256(3) + corner3(3) +', &
!                    'PERM(origin_mod_256(4) + corner3(4)) > 512'
!            print*, 'origin_mod_256(2) + corner3(2) + PERM(origin_mod_256(3) + corner3(3) +', &
!                    'PERM(origin_mod_256(4) + corner3(4))',&
!                    origin_mod_256(2) + corner3(2) + PERM(origin_mod_256(3) + corner3(3) + &
!                    PERM(origin_mod_256(4) + corner3(4)))
!            print*, 'origin_mod_256(2)', origin_mod_256(2)
!            print*, 'corner3(2)', corner3(2)
!            print*, 'PERM(origin_mod_256(3) + corner3(3) + PERM(origin_mod_256(4) + corner3(4))',&
!                    PERM(origin_mod_256(3) + corner3(3) + PERM(origin_mod_256(4) + corner3(4)))
!            print*, 'origin_mod_256(3)', origin_mod_256(3)
!            print*, 'corner3(3)', corner3(3)
!            print*, 'PERM(origin_mod_256(4) + corner3(4))', PERM(origin_mod_256(4) + corner3(4))
!            print*, 'origin_mod_256(4)', origin_mod_256(4)
!            print*, 'corner3(4)', corner3(4)
!        end if
!        if (origin_mod_256(2) + corner3(2) + PERM(origin_mod_256(3) + corner3(3) + &
!            PERM(origin_mod_256(4) + corner3(4))) < 1) then
!            print*, 'origin_mod_256(2) + corner3(2) + PERM(origin_mod_256(3) + corner3(3) +', &
!                    'PERM(origin_mod_256(4) + corner3(4)) < 1'
!            print*, 'origin_mod_256(2) + corner3(2) + PERM(origin_mod_256(3) + corner3(3) +', &
!                    'PERM(origin_mod_256(4) + corner3(4))',&
!                    origin_mod_256(2) + corner3(2) + PERM(origin_mod_256(3) + corner3(3) + &
!                    PERM(origin_mod_256(4) + corner3(4)))
!            print*, 'origin_mod_256(2)', origin_mod_256(2)
!            print*, 'corner3(2)', corner3(2)
!            print*, 'PERM(origin_mod_256(3) + corner3(3) + PERM(origin_mod_256(4) + corner3(4))',&
!                    PERM(origin_mod_256(3) + corner3(3) + PERM(origin_mod_256(4) + corner3(4)))
!            print*, 'origin_mod_256(3)', origin_mod_256(3)
!            print*, 'corner3(3)', corner3(3)
!            print*, 'PERM(origin_mod_256(4) + corner3(4))', PERM(origin_mod_256(4) + corner3(4))
!            print*, 'origin_mod_256(4)', origin_mod_256(4)
!            print*, 'corner3(4)', corner3(4)
!        end if
        gi1 = mod(PERM(origin_mod_256(1) + corner1(1) + &
                     PERM(origin_mod_256(2) + corner1(2) + &
                     PERM(origin_mod_256(3) + corner1(3) + &
                     PERM(origin_mod_256(4) + corner1(4))))), 32)
        gi2 = mod(PERM(origin_mod_256(1) + corner2(1) + &
                        PERM(origin_mod_256(2) + corner2(2) + &
                        PERM(origin_mod_256(3) + corner2(3) + &
                        PERM(origin_mod_256(4) + corner2(4))))), 32)
        gi3 = mod(PERM(origin_mod_256(1) + corner3(1) + &
                        PERM(origin_mod_256(2) + corner3(2) + &
                        PERM(origin_mod_256(3) + corner3(3) + &
                        PERM(origin_mod_256(4) + corner3(4))))), 32)
        gi4 = mod(PERM(origin_mod_256(1) + 1 + &
                        PERM(origin_mod_256(2) + 1 + &
                        PERM(origin_mod_256(3) + 1 + &
                        PERM(origin_mod_256(4) + 1)))), 32)
        t0 = 0.6 - sum(v0**2)
        if (t0 < 0.0) then
            n0 = 0.0
        else
            n0 = t0**4 * dot_product(v0, GRAD4(gi0+1, :))
        end if
        t1 = 0.6 - sum(v1**2)
        if (t1 < 0.0) then
            n1 = 0.0
        else
            n1 = t1**4 * dot_product(v1, GRAD4(gi1+1, :))
        end if
        t2 = 0.6 - sum(v2**2)
        if (t2 < 0.0) then
            n2 = 0.0
        else
            n2 = t2**4 * dot_product(v2, GRAD4(gi2+1, :))
        end if
        t3 = 0.6 - sum(v3**2)
        if (t3 < 0.0) then
            n3 = 0.0
        else
            n3 = t3**4 * dot_product(v3, GRAD4(gi3+1, :))
        end if
        t4 = 0.6 - sum(v4**2)
        if (t4 < 0.0) then
            n4 = 0.0
        else
            n4 = t4**4 * dot_product(v4, GRAD4(gi4+1, :))
        end if
        noise = 27.0 * (n0 + n1 + n2 + n3 + n4)
    end function snoise4d
end module simplex_noise_mod
