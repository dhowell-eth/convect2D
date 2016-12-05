% Remember full file needs to be transposed once subsetted
full = load('5/T_field.dat');
full_phi = load('5/phi_field.dat');
full_w = load('5/w_field.dat');

times = unique(full(:,1));

data_times = full(:,1);
full_data = full(:,2:end);
full_data_phi = full_phi(:,2:end);
full_data_w = full_w(:,2:end);

figure(1)
h1 = subplot(3,1,1);
h2 = subplot(3,1,2);
h3 = subplot(3,1,3);
%% 
for i=1:length(times)
   this_time = times(i);
   test = find(data_times==times(i));
   this_data = full_data(test,:)';
   this_data_phi = full_data_phi(test,:)';
   this_data_w = full_data_w(test,:)';
   
   contourf(h1,this_data);
   title(h1,strcat('Temperature: t= ',num2str(this_time)));
   contourf(h2,this_data_w);
   title(h2,strcat('Omega: t= ',num2str(this_time)));
   contourf(h3,this_data_phi);
   title(h3,strcat('Stream Function: t= ',num2str(this_time)));
   pause(1.5);
end

