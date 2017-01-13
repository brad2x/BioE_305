%%Fraction of Correct PCR Products

cycles = 30; % how many cycles of PCR to simulate?
construct_length = 1000;
error_rates = [4.4 36.7 330] * 10^(-7) * construct_length; % TODO: replace ? with a number
% such that this vector holds the error rates per DNA molecule per PCR cycle
initial_dsDNA = 10^10;

num_polmerases = length(error_rates);
cycle_num = 0:1:cycles; % start with cycle 0, then cycle 1, etc
frac_correct = ones(num_polmerases, cycles + 1); % a matrix where each row is a
% polymerase and each column is the fraction of correct products at a cycle

correct_dsDNA = zeros(num_polmerases, cycles); % a matrix where each row is a
% polymerase and each column is the number of correct dsDNA at a cycle
for i = 1:num_polmerases
    correct_dsDNA(i, 1) = initial_dsDNA;
end
incorrect_dsDNA = zeros(num_polmerases, cycles); % dsDNA with both strands being incorrect
half_correct_dsDNA = zeros(num_polmerases, cycles); % dsDNA with exactly one strand being incorrect
total_dsDNA = zeros(num_polmerases, cycles);

correct_ssDNA = zeros(num_polmerases, cycles);
incorrect_ssDNA = zeros(num_polmerases,cycles);
total_ssDNA = zeros(num_polmerases,cycles);

for j = 2:(cycles + 1) % iterate over cycles
    % At each cycle, dsDNA is melted into ssDNA and each undergoes replication.
    % A correct ssDNA might be the basis of making an incorrect ssDNA.
    % An incorrect ssDNA will be the basis of making an incorrect ssDNA.
    % Remember that we always have any ssDNA from the previous cycle.
    correct_ssDNA(:,j) = 2*correct_dsDNA(:,j-1) + half_correct_dsDNA(:,j-1);
    incorrect_ssDNA(:,j) = 2*incorrect_dsDNA(:,j-1) + half_correct_dsDNA(:,j-1);
    total_ssDNA(:,j) = correct_ssDNA(:,j) + incorrect_ssDNA(:,j);
    
    for i = 1:num_polmerases % iterate over polymerases
        % TODO: update the correct_dsDNA, incorrect_dsDNA, half_correct_dsDNA,
        % frac_correct matrices
        
        half_correct_dsDNA(i,j) = (error_rates(i))*correct_ssDNA(i,j);
        correct_dsDNA(i,j) = (1 - error_rates(i))*correct_ssDNA(i,j);
        incorrect_dsDNA(i,j) = incorrect_ssDNA(i,j);
        total_dsDNA(i,j) = correct_dsDNA(i,j)+half_correct_dsDNA(i,j)+incorrect_dsDNA(i,j);
        
        frac_correct(i,j) = (correct_dsDNA(i,j))/(half_correct_dsDNA(i,j) + incorrect_dsDNA(i,j) + correct_dsDNA(i,j));
    end
end

plot(cycle_num, frac_correct, 'o-') % plot the results
xlabel('Cycle Number'); ylabel('Fraction of Correct PCR Products');
legend('Phusion','Pfu','Taq'); 
title('Fraction of Correct PCR Products');
% TODO: label axes
axis([0 cycles 0 1])

%%