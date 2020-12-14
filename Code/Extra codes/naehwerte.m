function [Kalorien, Protein, Kohlenhydrate, Zucker, Fett, Ballaststoffe, Salz] = naehwerte(inputpatient)

%fddb.com information of food
% this function stores the information with the 
%nutritional values of the products

dbstop if error


folder= '/media/raid/Alex/Stimuli/Used Stimuli/';

folder=  [ folder , inputpatient,'sat1/'];
patdir=  [ folder];
folder=  [ folder, inputpatient,'sat1/'];




listdir = dir([folder, '/*.jpg']);
if numel(listdir) ~= 20
    error('not 20 stimuli in folder')
    
end

% check if strictly increasing the name
clear ff
for ll=1:numel(listdir)
    ff(ll,:) = (listdir(ll).name(1:3));
end

stimID =  str2num(ff);




if all(diff(stimID)>0) % check if strictly increasing the name; otherwise error
    
    for stim = 1:numel(listdir)
        
        switch stimID(stim)
            
            
            case 1
                
                %% 1
                %Balisto
                
                
                Kalorien(stim,1) = 503;
                Protein(stim,1) = 7;
                Kohlenhydrate(stim,1) = 61;
                Zucker(stim,1) = 43.8;
                Fett(stim,1) = 24.9;
                Ballaststoffe(stim,1) = 3.2;
                % Broteinheiten = 5.1
                % Wassergehalt = 16
                Salz(stim,1) = 0.0686;
                
                
                
            case 2
                %% 2
                
                %Bounty
                
                
                Kalorien(stim,1) = 488;
                Protein(stim,1) = 3.7 ;
                Kohlenhydrate(stim,1) = 58.9 ;
                Zucker(stim,1) = 48.2 ;
                Fett(stim,1) = 25.7 ;
                Ballaststoffe(stim,1) = 1.6 ;
                Broteinheiten = 4.9;
                Wassergehalt = 10;
                Salz(stim,1) = 0.25 ;
                
                
                
            case 4
                %% 4
                
                %corny
                
                Kalorien(stim,1) = 432;
                Protein(stim,1) = 5.5;
                Kohlenhydrate(stim,1) = 68.3;
                Zucker(stim,1) = 35.3;
                Fett(stim,1) = 14.2 ;
                Ballaststoffe(stim,1)  = 0 ;
                Broteinheiten = 5.7;
                Salz(stim,1) = 0.51 ;
                
                
            case 8
                %% 8
                %  kinder riegel
                
                
                
                Kalorien(stim,1) = 566 ;
                Protein(stim,1) = 8.7 ;
                Kohlenhydrate(stim,1) = 53.5 ;
                Zucker(stim,1) = 53.3 ;
                Fett(stim,1) = 35;
                Ballaststoffe(stim,1) = 0.9 ;
                Broteinheiten = 4.5;
                Salz(stim,1) = 0.313 ;
                
                
                E = 1.5 ;
                B2 = 0.42 ;
                B12 = 0.72 ;
                
                
                Magnesium = 45 ;
                Kalium = 400 ;
                Kalzium = 309;
                Phosphor = 245;
                
                
                
            case 09
                %% 9
                
                % kit kat chunky
                
                
                Kalorien(stim,1) = 527 ;
                Protein(stim,1) = 6.7 ;
                Kohlenhydrate(stim,1) = 61.1 ;
                Zucker(stim,1) = 51.3 ;
                Fett(stim,1) = 28.2 ;
                Ballaststoffe(stim,1) = 1.3 ;
                Broteinheiten = 5.1;
                Salz(stim,1) = 0.2 ;
                
                
            case 10
                %% 10
                
                %Leibniz pick up
                
                
                Kalorien(stim,1) = 516 ;
                Protein(stim,1) = 6.3 ;
                Kohlenhydrate(stim,1) = 61 ;
                Zucker(stim,1) = 34 ;
                Fett(stim,1) = 27 ;
                Ballaststoffe(stim,1) = 2.1 ;
                Broteinheiten = 5.1;
                Salz(stim,1)  = 0.5 ;
                
                
            case 12
                %% 12
                % mars
                
                
                
                Kalorien(stim,1) = 449 ;
                Protein(stim,1)  = 4 ;
                Kohlenhydrate(stim,1) = 70.2 ;
                Zucker(stim,1) = 61.7 ;
                Fett(stim,1) = 16.6 ;
                Ballaststoffe(stim,1) = 1.2 ;
                Broteinheiten = 5.9;
                Wassergehalt = 20;
                Salz(stim,1)  = 0.39 ;
                
                
                
            case 18
                %% 18
                % Role
                
                
                Kalorien(stim,1) = 480 ;
                Protein(stim,1) = 4.2 ;
                Kohlenhydrate(stim,1) = 68.3 ;
                Zucker(stim,1) = 59 ;
                Fett(stim,1) = 20.8 ;
                Ballaststoffe(stim,1)  = 1.1 ;
                Broteinheiten = 5.7;
                Salz(stim,1)= 0.3 ;
                
                
            case 19
                %% 19
                
                % snickers
                Kalorien(stim,1) = 481 ;
                Protein(stim,1) = 8.6 ;
                Kohlenhydrate(stim,1) = 60.5 ;
                Zucker(stim,1) = 51.8 ;
                Fett(stim,1) = 22.5 ;
                Ballaststoffe(stim,1) = 1.1 ;
                Broteinheiten = 5;
                Wassergehalt = 5;
                Salz(stim,1) = 0.63 ;
                
            case 20
                %% 20
                % twix
                
                
                Kalorien(stim,1) = 494;
                Protein(stim,1) = 4.5 ;
                Kohlenhydrate(stim,1) = 64.6 ;
                Zucker(stim,1) = 48.6 ;
                Fett(stim,1) = 23.9 ;
                Ballaststoffe(stim,1) = 1.5 ;
                Broteinheiten = 5.4;
                Wassergehalt = 20;
                Salz(stim,1) = 0.41 ;
                
                
                
            case 101
                
                %% 101
                
                % bifi carazza
                
                
                Kalorien(stim,1) = 393;
                Protein(stim,1) = 12 ;
                Kohlenhydrate(stim,1) =  38 ;
                Zucker(stim,1) = 9.3;
                Fett(stim,1) = 22;
                Ballaststoffe(stim,1) = 2 ;
                Broteinheiten = 3.2;
                Wassergehalt = 30;
                Salz(stim,1) = 2.032 ;
                
                
                
                
            case 102
                %% 102
                
                % bifi mussmit
                
                Kalorien(stim,1) = 510 ;
                Protein(stim,1) = 25 ;
                Kohlenhydrate(stim,1) = 0.9;
                Zucker(stim,1) = 0.9 ;
                Fett(stim,1) = 45 ;
                Ballaststoffe(stim,1) = 0.1 ;
                Broteinheiten = 0.1;
                Wassergehalt = 15;
                
                Salz(stim,1) = 4.0132 ;
                
                
            case 104
                % 104
                
                % crocantelli alle olive
                
                Kalorien(stim,1) = 484;
                Protein(stim,1) = 11.7;
                Kohlenhydrate(stim,1) = 74.2;
                Zucker(stim,1) =  1.7 ;
                Fett(stim,1) = 14.4;
                Ballaststoffe(stim,1) = 5.0 ;
                %                 Broteinheiten =;
                %                 Wassergehalt = ;
                Salz(stim,1) = 2.5;
                
                
                
            case 107
                
                %% 107
                % Zwiebliringe
                
                Kalorien(stim,1) = 515;
                Protein(stim,1)	 = 7 ;
                Kohlenhydrate(stim,1) =	57;
                Zucker(stim,1) = 3.6 ;
                Fett(stim,1) = 28 ;
                Ballaststoffe(stim,1) = 3.9 ;
                
                Salz(stim,1) = 1.5 ;
                
                
                
                
            case 109
                %% 109
                %glod fischli sesam
                
                Kalorien(stim,1) = 483 ;
                Protein(stim,1) = 11 ;
                Kohlenhydrate(stim,1) = 49 ;
                Zucker(stim,1) = 6.2 ;
                Fett(stim,1) = 26 ;
                Ballaststoffe(stim,1) = 4.6 ;
                Broteinheiten = 4.1;
                Salz(stim,1) = 2.54 ;
                
                
                
                
            case 111
                %% 111
                % lays light gesalzen
                
                Kalorien(stim,1)  = 490 ;
                Protein(stim,1) = 7 ;
                Kohlenhydrate(stim,1) = 64 ;
                Zucker(stim,1)	 = 0.3 ;
                Fett(stim,1) = 22 ;
                Ballaststoffe(stim,1) = 4.5 ;
                Broteinheiten = 5.3;
                Salz(stim,1) = 1.43 ;
                
                
                
            case 113
                %% 113
                % erdnuß locken
                
                Kalorien(stim,1) = 500 ;
                Protein(stim,1) = 13 ;
                Kohlenhydrate(stim,1) = 56 ;
                Zucker(stim,1) = 2.3 ;
                Fett(stim,1) = 24 ;
                Ballaststoffe(stim,1) = 4.1 ;
                Broteinheiten = 4.7;
                Wassergehalt = 5;
                Salz(stim,1) = 1.9304 ;
                
                
                
            case 115
                %% 115
                % Saltletts
                
                Kalorien(stim,1) = 399 ;
                Protein(stim,1) = 12 ;
                Kohlenhydrate(stim,1) = 69 ;
                Zucker(stim,1) = 2.7 ;
                Fett(stim,1) = 7.5 ;
                Ballaststoffe(stim,1) = 3.7;
                Broteinheiten = 5.8;
                Salz(stim,1) = 4.1;
                
                
                
            case 117
                %%  117
                % tuc
                
                
                Kalorien(stim,1) = 479;
                Protein(stim,1) = 9.2;
                Kohlenhydrate(stim,1) = 64;
                Zucker(stim,1) = 7.8;
                Fett(stim,1) =  20;
                Ballaststoffe(stim,1) = 2.9;
                
                Salz(stim,1) = 2.7;
                
                
                
            case 119
                %% 119
                % ültje erdnüsse
                
                Kalorien(stim,1) = 615 ;
                Protein(stim,1) = 25 ;
                Kohlenhydrate(stim,1) = 14 ;
                
                Zucker(stim,1) = 5.4 ;
                
                Fett(stim,1) = 50 ;
                Ballaststoffe(stim,1) = 6.7 ;
                Broteinheiten = 1.2;
                Salz(stim,1) = 1.016 ;
                
                
                
            case 120
                %% 120
                % tortillA CHIPS CHILI
                
                
                Kalorien(stim,1) = 483 ;
                Protein(stim,1) = 5.9 ;
                Kohlenhydrate(stim,1) = 63 ;
                Zucker(stim,1)  = 2.9 ;
                Fett(stim,1)  = 22 ;
                Ballaststoffe(stim,1) = 4.3 ;
                Broteinheiten = 5.3;
                Salz(stim,1) = 2 ;
                
                
                
        end
        
    end
end
