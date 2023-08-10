% subroutine SwashBCboundwave ( bcfour, nfreq, xp, yp, ibgrpt, swd, wdir, rsgn, vdir, shape, ibloc )
% !
% !   --|-----------------------------------------------------------|--
% !     | Delft University of Technology                            |
% !     | Faculty of Civil Engineering                              |
% !     | Environmental Fluid Mechanics Section                     |
% !     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
% !     |                                                           |
% !     | Programmers: The SWASH team                               |
% !   --|-----------------------------------------------------------|--
% !
% !
% !     SWASH (Simulating WAves till SHore); a non-hydrostatic wave-flow model
% !     Copyright (C) 2010-2018  Delft University of Technology
% !
% !     This program is free software; you can redistribute it and/or
% !     modify it under the terms of the GNU General Public License as
% !     published by the Free Software Foundation; either version 2 of
% !     the License, or (at your option) any later version.
% !
% !     This program is distributed in the hope that it will be useful,
% !     but WITHOUT ANY WARRANTY; without even the implied warranty of
% !     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% !     GNU General Public License for more details.
% !
% !     A copy of the GNU General Public License is available at
% !     http://www.gnu.org/copyleft/gpl.html#SEC3
% !     or by writing to the Free Software Foundation, Inc.,
% !     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
% !
% !
% !   Authors
% !
% !    1.00: Dirk Rijnsdorp
% !
% !   Updates
% !
% !    1.00, October 2012: New subroutine
% !
% !   Purpose
% !
% !   Computes second order bound long wave to be added to Fourier series along open boundaries
% !
% !   Method
% !
% !   Details on the computation of the interaction coefficient of two primary wave components can be found in
% !
% !   K. Hasselmann
% !   On the non-linear energy transfer in a gravity-wave spectrum
% !   J. Fluid Mech., vol 12, 481-500, 1962
% !
% !   Modules used
% !
%     use ocpcomm4
%     use SwashCommdata3
%     use m_bndspec
% !
%     implicit none
% !
% !   Argument variables
% !
%     integer, intent(in)         :: ibgrpt ! actual boundary grid point
%     integer, intent(in)         :: ibloc  ! actual counter for boundary point
%     integer, intent(in)         :: nfreq  ! number of frequencies
%     integer, intent(in)         :: shape  ! spectral shape
%                                           ! = 1; Pierson Moskowitz
%                                           ! = 2; Jonswap
%                                           ! = 3; TMA
%     !
%     real, intent(in)            :: rsgn   ! sign for indicating in- and outflowing depending on boundary
%                                           ! =+1; refers to inflowing at left and lower boundaries
%                                           ! =-1; refers to outflowing at right and upper boundaries
%     real, intent(in)            :: swd    ! still water depth
%     real, intent(in)            :: wdir   ! incident or peak wave direction with respect to problem coordinates
%     real, intent(in)            :: xp     ! x-coordinate of grid point
%     real, intent(in)            :: yp     ! y-coordinate of grid point
%     !
%     logical, intent(in)         :: vdir   ! indicates direction of in- or outcoming velocity on boundary
%                                           ! =.true.; u-velocity
%                                           ! =.false.; v-velocity
%     !
%     type(bfsdat), intent(inout) :: bcfour ! list containing parameters for Fourier series
% !
% !   Parameter variables
% !
%     real, parameter :: khmin = 0.4 ! minimum dimensionless depth
%     real, parameter :: khmax = 10. ! maximum dimensionless depth
% !
% !   Local variables
% !
%     integer, save :: ient  = 0    ! number of entries in this subroutine
%     integer, save :: inext = 0    ! wait for next boundary point to allocate
%     integer       :: j            ! loop counter
%     integer       :: k            ! loop counter
%     integer       :: l            ! frequency index of bound long wave
%     !
%     real          :: ampl1        ! amplitude of first primary wave component
%     real          :: ampl2        ! amplitude of second primary wave component
%     real          :: ampl3        ! amplitude of bound long wave component
%     real          :: cg3          ! group velocity of bound long wave
%     real          :: cosdt        ! cosine of difference angles of primary waves
%     real          :: df           ! increment in frequency space
%     real          :: Dp           ! interaction coefficient for second order potential amplitude as given by Eq. (4.3) of Hasselmann (1962)
%     real          :: Dz           ! interaction coefficient for second order surface elevation amplitude as given by Eq. (4.2) of Hasselmann (1962)
%     real          :: f1           ! frequency of first primary wave component
%     real          :: f2           ! frequency of second primary wave component
%     real          :: f3           ! frequency of bound long wave component
%     real          :: kwav1        ! wave number of first primary wave component
%     real          :: kwav2        ! wave number of second primary wave component
%     real          :: kwav3        ! wave number of bound long wave component
%     real          :: n            ! ratio of group and phase velocity
%     real          :: omega1       ! angular frequency of first primary wave component
%     real          :: omega2       ! angular frequency of second primary wave component
%     real          :: phase1       ! phase of first primary wave component
%     real          :: phase2       ! phase of second primary wave component
%     real          :: phase3       ! phase of bound long wave component
%     real          :: rval         ! auxiliary real
%     real          :: s            ! sign
%     real          :: T1           ! first term of right-hand side of Eq. (4.2) of Hasselmann (1962)
%     real          :: T2           ! second term of right-hand side of Eq. (4.2) of Hasselmann (1962)
%     real          :: theta1       ! wave direction of first primary wave component
%     real          :: theta2       ! wave direction of second primary wave component
%     real          :: theta3       ! wave direction of bound long wave component
% !
% !   Structure
% !
% !   Description of the pseudo code
% !
% !   Source text
% !
%     if (ltrace) call strace (ient,'SwashBCboundwave')
%     !
%     ! direction of each wave component in line with incident direction
%     ! in order to preserve symmetry at boundaries
%     !
%     if ( vdir ) then
%        if ( .not. sin(wdir) < 0. ) then
%           s = +1.
%        else
%           s = -1.
%        endif
%     else
%        if ( .not. cos(wdir) < 0. ) then
%           s = +1.
%        else
%           s = -1.
%        endif
%     endif
%     !
%     if ( inext < ibloc ) then
%        !
%        allocate(bcfour%zetab (nbgrpt,nfreq))
%        allocate(bcfour%fluxbu(nbgrpt,nfreq))
%        allocate(bcfour%fluxbv(nbgrpt,nfreq))
%        !
%        bcfour.zetab  = (0.,0.)
%        bcfour.fluxbu = (0.,0.)
%        bcfour.fluxbv = (0.,0.)
%        !
%        inext = ibloc
%        !
%     endif
%     !
%     df = ( bcfour%omega(2) - bcfour%omega(1) ) / pi2
%     !

pi2 = 2*pi;
%     floop: do j = 1, nfreq
for j=1:nfreq %floop
    % loop through all frequencies for j = 1:nfreq
    
    % get first primary wave component
    
    ampl1  = bcfour.ampl (j);
    omega1 = bcfour.omega(j);
    phase1 = bcfour.phase(j);
    
    theta1 = wdir + s*bcfour.theta(j);
    
    % calculate frequency of this first component
    
    f1 = omega1 / pi2;
    
    % calculate wave number of this component
    
    kwav1 = getk(f1,swd);
    
    %        call disprel ( swd, omega1, kwav1, rval, n )
    
    %     % % correct amplitude in case of TMA spectrum for shallow water
    %      if ( shape == 3 )
    %          ampl1 = ampl1 * omega1 * omega1 / ( grav * kwav1 * sqrt(2.*n) );
    %      end
    %
    % % in case of periodicity, wave direction must be corrected so that wave number is an integer multiple of 2pi/length with length the periodicity length
    %
    %        if ( bcperx ) then
    %           !
    %           rval = nint( kwav1*cos(theta1) / ( pi2/xclen ) ) * pi2/xclen / kwav1
    %           if ( rval > 1. ) then
    %              theta1 = acos ( rval - pi2/xclen / kwav1 )
    %           else if ( rval < -1. ) then
    %              theta1 = acos ( rval + pi2/xclen / kwav1 )
    %           else
    %              theta1 = acos ( rval )
    %           endif
    %           !
    %        else if ( bcpery ) then
    %           !
    %           rval = nint( kwav1*sin(theta1) / ( pi2/yclen ) ) * pi2/yclen / kwav1
    %           if ( rval > 1. ) then
    %              theta1 = asin ( rval - pi2/yclen / kwav1 )
    %           else if ( rval < -1. ) then
    %              theta1 = asin ( rval + pi2/yclen / kwav1 )
    %           else
    %              theta1 = asin ( rval )
    %           endif
    %           !
    %        endif
    
    % check this direction with respect to the normal of boundary
    % (must be within -80 degrees to 80 degrees)
    
    if  vdir
        if rsgn == 1
            if ( cos(theta1) <  0.174 )
                continue
            elseif ( cos(theta1) > -0.174 )
                continue
            end
        end
    else
        if rsgn == 1
            if ( sin(theta1) <  0.174 )
                continue
            elseif ( sin(theta1) > -0.174 )
                continue
            end
        end
    end
    
    if ( kwav1*swd > khmin .and. kwav1*swd < khmax .and. ampl1 ~= 0 )
        
        for k = j+1: nfreq
            
            % get second primary wave component
            
            ampl2  = bcfour.ampl (k);
            omega2 = bcfour.omega(k);
            phase2 = bcfour.phase(k);
            
            theta2 = wdir + s*bcfour.theta(k);
            
            % calculate frequency of this second component
            
            f2 = omega2 / pi2;
            
            % calculate wave number of this component
            
            %              call disprel ( swd, omega2, kwav2, rval, n )
            k2 = getk(f2,swd);
            
            %             % correct amplitude in case of TMA spectrum for shallow water
            %
            %             if ( shape == 3 )
            %                 ampl2 = ampl2 * omega2 * omega2 / ( grav * kwav2 * sqrt(2.*n) );
            %             end
            %
            %             % in case of periodicity, wave direction must be corrected so that wave number is an integer multiple of 2pi/length with length the periodicity length
            %
            %             if ( bcperx )
            %
            %                 rval = nint( kwav2*cos(theta2) / ( pi2/xclen ) ) * pi2/xclen / kwav2;
            %                 if ( rval > 1. )
            %                     theta2 = acos ( rval - pi2/xclen / kwav2 );
            %                 elseif ( rval < -1. )
            %                     theta2 = acos ( rval + pi2/xclen / kwav2 );
            %                 else
            %                     theta2 = acos ( rval );
            %                 end
            %
            %             elseif ( bcpery )
            %
            %                 rval = nint( kwav2*sin(theta2) / ( pi2/yclen ) ) * pi2/yclen / kwav2;
            %                 if ( rval > 1. )
            %                     theta2 = asin ( rval - pi2/yclen / kwav2 );
            %                 elseif ( rval < -1. )
            %                     theta2 = asin ( rval + pi2/yclen / kwav2 );
            %                 else
            %                     theta2 = asin ( rval );
            %                 end
            %
            %             end
            
            % check this direction with respect to the normal of boundary
            % (must be within -80 degrees to 80 degrees)
            
            if ( vdir )
                if ( rsgn == 1. )
                    if ( cos(theta2) <  0.174 )
                        continue
                        
                    else
                        if ( cos(theta2) > -0.174 )
                            continue
                        end
                    end
                    
                elseif ( rsgn == 1. )
                    if ( sin(theta2) <  0.174 )
                        continue
                    else
                        if ( sin(theta2) > -0.174 )
                            continue
                        end
                    end
                end
            end
            
            if ( kwav2*swd > khmin .and. kwav2*swd < khmax .and. ampl2 ~= 0. )
                
                % calculate frequency of bound long wave
                f3 = f2 - f1;
                % calculate direction of bound long wave 
                theta3 = atan( ( kwav2*sin(theta2) - kwav1*sin(theta1) ) / ( kwav2*cos(theta2) - kwav1*cos(theta1) ) );
                % calculate cosine of difference angles
                cosdt = cos( theta1 - theta2 + pi );
                % calculate bound wave number
                kwav3 = real(sqrt( dble(kwav1)^.2 + dble(kwav2).^2 + 2.*doubble(kwav1)*double(kwav2)*cosdt ));
                % calculate bound wave group velocity
                cg3 = pi2 * f3 / kwav3;
                % calculate bound wave phase
                phase3 = phase1 - phase2 + pi;
                % include phase shift related to wave direction and wave number
                phase3 = phase3 + kwav3 * ( cos(theta3)*xp + sin(theta3)*yp );
                % for computing the interaction coefficient for bound long wave, change sign of omega2
                omega2 = -omega2;
                % compute interaction coefficient as given by Eq. (4.3) of Hasselmann (1962)
                Dp = (omega1+omega2) * ( (omega1*omega2).^2 /(grav^2) - kwav1*kwav2*cosdt ) - 0.5 * ( omega1*kwav2.^2/(cosh(kwav2*swd).^2) + omega2*kwav1.^2/(cosh(kwav1*swd).^2) );
                % compute the time derivative of second order potential
                % (derived from Eq. (4.7) of Hasselmann (1962) and use relation between potential and surface elevation for transformation, as given next to Eq. (1.26) at page 485)
                T1 = Dp*grav*(omega1+omega2)/(omega1*omega2) / ( grav*kwav3*tanh(kwav3*swd) - (omega1+omega2).^2 );
                % compute interaction coefficient as given by Eq. (4.4) of Hasselmann (1962)
                % (use relation between potential and surface elevation for transformation, as given next to Eq. (1.26) at page 485)
                T2 = -grav*kwav1*kwav2 / (2.*omega1*omega2) * cosdt + 0.5 * ( omega1^2 + omega1*omega2 + omega2^2 ) / grav;
                % compute interaction coefficient for second order surface elevation amplitude, Eq. (4.2) of Hasselmann (1962)
                Dz = T1 + T2;
                % compute amplitude of bound wave
                ampl3 = abs(Dz) * ampl1 * ampl2;
                l = round(f3/df);
                % store surface elevation and mass flux of bound long wave
                bcfour.zetab(ibgrpt,l) = bcfour.zetab(ibgrpt,l) + ampl3 * exp( (0.,1.) * phase3 );
                
                if ( oned )
                    bcfour.fluxbu(ibgrpt,l) = bcfour.fluxbu(ibgrpt,l) + cg3 * ampl3 * exp( (0.,1.) * phase3 );
                else
                    bcfour.fluxbu(ibgrpt,l) = bcfour.fluxbu(ibgrpt,l) + cg3 * ampl3 * cos(theta3) * exp( (0.,1.) * phase3 );
                    bcfour.fluxbv(ibgrpt,l) = bcfour.fluxbv(ibgrpt,l) + cg3 * ampl3 * sin(theta3) * exp( (0.,1.) * phase3 );
                end
                
            end % kwav2
            
        end %sloop
        
    end %kwav1
    
end % floop

% end subroutine SwashBCboundwave
