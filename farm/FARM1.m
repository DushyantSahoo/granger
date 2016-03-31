arg_list = argv ();

fMRI_filename  = char(arg_list{1});
addpath("/usr/lib/granger/");

% Main function
	% Intialize 
    [intial_location,hdr1]=cbiReadNifti(fMRI_filename);

    y = intial_location;

    Size = size(y);
    a = Size(1);
    b = Size(2);
    c = Size(3);
    d = Size(4);
    itemp = 1;
	% Threshold filtered_func_data.nii.gz file and convert 3D array to 1D
    for i = 1:a
        for j = 1:b
            for k = 1:c
                if (mean(y(i,j,k,:)) > 6000)
                    number(:,itemp) = y(i,j,k,:);
                    itemp = itemp+1;
                end
            end
        end
    end
	% Intialization of paramters
    number = double(number);
    Size1 = size(number);
    a1 = Size1(1);
    b1 = Size1(2);
    Beta = zeros(b1,b1);
    beta3 = zeros(b1,1);
    beta4 = zeros(b1,1);
    beta5 = zeros(b1,1);
    beta6 = zeros(b1,1);
	% Granger Causality Analysis
    for i1 = 1:b1
        

        number1 = number;
        number1(:,i1) = [];
        temp1 = number(2:d,i1);
        temp2 = number1(1:(d-1),:);
        new = zscore(temp1) ;
        new = new/norm(new,2);
        new1 = zscore(temp2);
        new1 = new1/norm(new1(:,2),2);
	% Lasso function
   
        [beta ,steps,G,residuals,error,drop] = lasso(new1, new, 0, false,false,10^(-9),str2double(arg_list{2}));

        if i1 == 1
            Beta(i1,1) = 0;
            Beta(i1,2:b1) = beta; 
        elseif i1 == b1
            Beta(i1,1:(b1-1)) = beta;
            Beta(i1,b1) = 0;
        else
            Beta(i1,1:(i1-1)) = beta(1:(i1-1));
            Beta(i1,i1) = 0;
            Beta(i1,(i1+1):b1) = beta((i1):(b1-1));
        end
        
         beta3(i1) = steps;
         beta4(i1) = 1-residuals;
        
        
    end
    beta1 = sum(abs(Beta),1); % prediction power
    beta2 = sum(abs(Beta),2); % effect of all voxels on a particular voxel
    invbeta1 = zeros(a,b,c);
    invbeta2 = zeros(a,b,c);
    invbeta3 = zeros(a,b,c);
    invbeta4 = zeros(a,b,c);
    itemp = 1;
	% Converting all 1D array to 3D array
    for i = 1:a
        for j = 1:b
            for k = 1:c
                if (mean(y(i,j,k,:)) > 6000)
                    invbeta1(i,j,k) = beta1(itemp);
                    invbeta2(i,j,k) = beta2(itemp);
                    invbeta3(i,j,k) = beta3(itemp);
                    invbeta4(i,j,k) = beta4(itemp);
                    itemp = itemp+1;
                end
            end
        end
    end
	% Save all the files
    csvwrite(strcat(char(arg_list{3}),'/beta.csv'),invbeta1);
    cbiWriteNifti(invbeta1,strcat(char(arg_list{3}),'/beta1.nii.gz'),'d',size(invbeta1));
    cbiWriteNifti(invbeta1,strcat(char(arg_list{3}),'/beta2.nii.gz'),'d',size(invbeta2));
    cbiWriteNifti(invbeta1,strcat(char(arg_list{3}),'/steps.nii.gz'),'d',size(invbeta3));
    cbiWriteNifti(invbeta1,strcat(char(arg_list{3}),'/residuals.nii.gz'),'d',size(invbeta2));
    clear
