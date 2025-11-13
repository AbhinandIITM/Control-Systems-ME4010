clc; clear ;

% Revised script structure: system definition first, observer prep next, and
% simulation/plotting at the end. Computation sequences remain unchanged.

m1=3; m2=5; c1=2; c2=5; k1=4; k2=12;
u = 0;
A = [0 1 0 0;
    -(k1+k2)/m1 -(c1+c2)/m1 k2/m1 c2/m1;
    0 0 0 1;
    k2/m2 c2/m2 -k2/m2 -c2/m2];
B = [0; 0; 0; 1/m2];
C = [1 0 0 0];

% Partitioned matrices for the minimum-order observer design
A_aa = [0];
A_ab = [1 0 0];
A_ba = [-(k2+k1)/m1 0 k2/m2]';
A_bb = [-(c1+c2)/m1 k2/m1 c2/m1;
     0 0 1;
     c2/m2 -k2/m2 -c2/m2];

B_a = [0];
B_b = [0; 0; 1/m2];

C_a = [1];
C_b = [0 0 0];

% Observability check for reduced subsystem
rank(obsv(A_bb,A_ab));

eig(A);

% Desired dynamics for reduced-order observer
p_o = [-5.9 -5.8 -1];
Ke = place(A_bb', A_ab', p_o)';

A_hat = A_bb - Ke*A_ab;
B_hat = A_hat*Ke + A_ba - Ke*A_aa;
F_hat = B_b - Ke*B_a;
C_hat = [zeros(1, 3);
         eye(3)];
D_hat = [1;
         Ke];

% Running the true system
t_sim = 0:0.01:10;
x0 = [1; 0; -1; 0];
[t, x_true] = ode45(@(t, x) smd(x, A, B, u), t_sim, x0);
y_true = (C * x_true')';

% Adding noise to the measurement
y_noise = awgn(y_true, 30, 'measured');

% Observer state initialization
z0 = [2; 1.5; 1];

% Integrating observer equations
[t_o, z] = ode45(@(t_o, z) moo_dyn(t_o, z, A_hat, B_hat, F_hat, u, y_true, t), t_sim, z0);
[t_o_noise, z_noise] = ode45(@(t_o, z) moo_dyn(t_o, z, A_hat, B_hat, F_hat, u, y_noise, t), t_sim, z0);

% State reconstruction
x_hat = (C_hat * z')' + (D_hat * y_true')';
x_hat_noise = (C_hat * z_noise')' + (D_hat * y_noise')';

% Plotting for Mass 1
figure('Name', 'Observer: Mass 1 States (x1, x2)', 'NumberTitle', 'off'); %[output:1701b86f]
subplot(2,2,1); %[output:1701b86f]
plot(t, x_true(:,1), 'k', t_o, x_hat(:,1), 'b', t_o_noise, x_hat_noise(:,1), 'r--','LineWidth',1.2); %[output:1701b86f]
xlabel('Time (s)'); ylabel('Displacement (m)'); title('Mass 1 Displacement'); grid on; %[output:1701b86f]
legend('True','Estimated','Noisy'); %[output:1701b86f]

subplot(2,2,2); %[output:1701b86f]
plot(t, x_true(:,2), 'k', t_o, x_hat(:,2), 'b', t_o_noise, x_hat_noise(:,2), 'r--','LineWidth',1.2); %[output:1701b86f]
xlabel('Time (s)'); ylabel('Velocity (m/s)'); title('Mass 1 Velocity'); grid on; %[output:1701b86f]
legend('True','Estimated','Noisy'); %[output:1701b86f]

subplot(2,2,3); %[output:1701b86f]
plot(t, x_true(:,1)-x_hat(:,1), 'b', t, x_true(:,1)-x_hat_noise(:,1), 'r--','LineWidth',1.2); %[output:1701b86f]
xlabel('Time (s)'); ylabel('Error'); title('x1 Estimation Error'); grid on; %[output:1701b86f]
legend('No noise','With noise'); %[output:1701b86f]

subplot(2,2,4); %[output:1701b86f]
plot(t, x_true(:,2)-x_hat(:,2), 'b', t, x_true(:,2)-x_hat_noise(:,2), 'r--','LineWidth',1.2); %[output:1701b86f]
xlabel('Time (s)'); ylabel('Error'); title('x2 Estimation Error'); grid on; %[output:1701b86f]
legend('No noise','With noise'); %[output:1701b86f]

sgtitle('Observer: Mass 1 Analysis'); %[output:1701b86f]

% Plotting for Mass 2
figure('Name','Observer: Mass 2 States (x3, x4)','NumberTitle','off'); %[output:352ed959]
subplot(2,2,1); %[output:352ed959]
plot(t, x_true(:,3), 'k', t_o, x_hat(:,3), 'b', t_o_noise, x_hat_noise(:,3), 'r--','LineWidth',1.2); %[output:352ed959]
xlabel('Time (s)'); ylabel('Displacement (m)'); title('Mass 2 Displacement'); grid on; %[output:352ed959]
legend('True','Estimated','Noisy'); %[output:352ed959]

subplot(2,2,2); %[output:352ed959]
plot(t, x_true(:,4), 'k', t_o, x_hat(:,4), 'b', t_o_noise, x_hat_noise(:,4), 'r--','LineWidth',1.2); %[output:352ed959]
xlabel('Time (s)'); ylabel('Velocity (m/s)'); title('Mass 2 Velocity'); grid on; %[output:352ed959]
legend('True','Estimated','Noisy'); %[output:352ed959]

subplot(2,2,3); %[output:352ed959]
plot(t, x_true(:,3)-x_hat(:,3), 'b', t, x_true(:,3)-x_hat_noise(:,3), 'r--','LineWidth',1.2); %[output:352ed959]
xlabel('Time (s)'); ylabel('Error'); title('x3 Estimation Error'); grid on; %[output:352ed959]
legend('No noise','With noise'); %[output:352ed959]

subplot(2,2,4); %[output:352ed959]
plot(t, x_true(:,4)-x_hat(:,4), 'b', t, x_true(:,4)-x_hat_noise(:,4), 'r--','LineWidth',1.2); %[output:352ed959]
xlabel('Time (s)'); ylabel('Error'); title('x4 Estimation Error'); grid on; %[output:352ed959]
legend('No noise','With noise'); %[output:352ed959]

sgtitle('Observer: Mass 2 Analysis'); %[output:352ed959]

function dx = smd(x, A, B, u)
    dx = A*x + B*u;
end

function dz = moo_dyn(t_o, z, A_hat, B_hat, F_hat, u, y_vec, t_vec)
    y_interp = interp1(t_vec, y_vec, t_o);
    dz = A_hat*z + B_hat*y_interp + F_hat*u;
end


%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":26.1}
%---
%[output:1701b86f]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAMUAAAB2CAYAAAB4ZlcrAAAAAXNSR0IArs4c6QAAHoxJREFUeF7tXQmQFMXSzllYWBX4FXTl1MUDvOWJD\/BAMBQQUcAbjVAgxBAVRUPAEzkUBETCGw01EENQvG9\/9XmgghAIavzhESqyHE+R83EI67I787+vdnPMqemjurenp2e3+8ULgamurvoyv8qsKzORSqVSFD8xAjECaQQSMSlibYgRyEQgJkWsETECGgIxKWKViBGISRHrQIyAMwKxpYg1JEYgSEsxffp0euKJJ9JVHn\/88fT0009Ty5Yt0\/+GMnhuueWWegv+ggUL6IUXXsjqO3d4y5YtdOWVV1KPHj3yisPPP\/9MU6dOpfvvvz9DRnaCQb9uu+02uvfee+mSSy4JVH5oy7hx42jGjBl0+OGH+64bbVy8eDFNmzaN9tprL9\/1yBd9WQoW8kEHHZTRmCVLltBll11G8+fPVwqApyGR4uijj6Y77rgjSziPPPIIffTRR3klBZRw+PDhdMABB9iSVyoGy7hFixa0fft2o3e8aGRQpPDyTdOyvkjhpOggxpgxY2jOnDlqBGgopMBo1bZtWzr\/\/PMzRj4o15NPPkk7duyg5s2b58VS8Ijfs2dPYwWHHCE7WHjIc+bMmemBzlS5nMrVK1K4dWb37t1066230sknn6xMLoCFQuD\/b731lsLp6quvTisHl+ffoFhMKJTVf5cuGrsDXbt2pVmzZtGAAQNU\/aeeemqGuWcBs2sn3T79e\/htn332USP7t99+q6zewoULHV1ANuGDBg2ijRs3Znwbv+EpLy\/PqEN3PWU73DBhi8xKJ\/HUFRGkvPnmm+n222+nzZs3K3noLq6V8vJgdsMNNyh5on3SBebf8S670Oeee26G5+DUR6lHr776aha+0i0CHnA\/IQ88Ugd094kHAO6TH9fPs6XQFcwKUNnQhx56SIHGLhWb5SFDhqRJgzoYcCswuCzKSf8dQoZLAGHI93X\/3knAVpYNBJXEdBsZuc3sI0+YMEH57BDmlClTaOjQoSQFrwuSSYDvwDcGZnaYrFu3LsMXZzzRf3ZZ7dprIju8y64WWwer91jh7eTq1kfZD52scmAdOHBgxiCL9knvQ34HpJGkdxvA7XDyRYrnn3\/ecWKjk0IKGH9mkOFrY6Klj0LcWJTTvyUBO+GEExQppGnXlUT+vVWrVpYTTQky\/vzbb795mrjp\/e3Vq5dSUAhl7ty5ap6hK7ouEK7jrrvuosmTJztiIt1TN8LK301JoSu0PpCxYuo4ubnKsl5JCsgFloCJLZWZf5MDo+yTrPPNN990XPAwxcoXKdxMsK4kZWVlGS6F7DRGCUzO8eiujG4KZadgFkEKqxUMKRypCD\/99FP6WzpA7IK4CdbNMmK0grsFAbPrxG6k1eDAfcdv7H6gDjtMdNfKyXXS22pCCt395Tr0wcIKJ6t\/01097qNu8exGf6woST3QVzid3Cw\/rhP665kUbibJak7hRApejpPCZnLA5XAate3aIoWPSS5\/323pVDfNpiOLLphJkyYpss6ePVu5TvqCAwtZDgJWS4tWmMjlS6lwJuQwIYWuxBID2V43Urj1UScFzw\/vuece5T3wnFR+ny0WBg0mx4cffmi5JCuJ5JUcnknhpjgmq09OwpGkwrec1v\/tSMHgwR+FSeW1cL1tVkpfV0uBkQ11lJaW0sqVK9NLtPq8Rhe6k9tmN3pLN9PNeku31WmibdcOfa\/FiRQ8OXfqo04K7mOXLl0yZGYlI+kSr1692nGfws8+hi9SeN2nsJto8yRKzims\/Em56SUngfA37TaAeKSQKyL6hBYK7DbxN7EWdisgcoRymuzLtlrNKXR3U59TmBLZzVK4TdrxHdQBUsEC6+6gaR+xmKCTAnVZycxqPiP7IS2FPqdwG0zsZOuLFFyZ7vPb7WjrS7JSWaRJ5Hrl5p\/T706uHJPn+uuvz5jP6D45vmmlvFbLj3a78joprNolFVfvE4iL5VwQAqtePLnkJUi0UWKi464vhdoJ240Ubu6lHJCslqm99BFt1Ac0O5nxv8OV1ueeOvb6MrCJW6njVSdSmIyicZkYAVME3OarpvXUtVxMiroiGL8fGAIY9bHJme9zcjEpAhNpXJFfBNidxPsmu+1+v2P6XkwKU6Ticg0GgZgUDUbUcUdNEYhJYYpUXK7BIBCTosGIOu6oKQIxKUyRiss1GARiUjQYUccdNUUgJoUpUnG5BoNATIoGI+q4o6YIxKQwRSou12AQsCWFPHSmR+fAqVf9QlAuEUNbHn74YcsrolYHBtEW\/ZAfn+6U4Xe8tFmeEPVbh5fvOZX1IxvccuzXr1+dwslYtclONnzwcs2aNRm71HYnrGXdVlFhvGKn14FzVe+\/\/z6NGjXKtSojUvBJQ7dLL65f81GATz3akdAqphID4vVyiVPzokoKE9k4DSo+RJJ+xU02TF45qJrIJghSyH55jbvlSgqERdl3333VnWWcgced6EMPPVRdoLE65iyP6tpdI9RHdzvllceAvZACgOhn\/6WlkPXysWu8g6gViHGEB0ej9ZtmdnXII\/NO10UlHvxdjlSBmFEc8QT1XXXVVelRTeKDOnBgrqioiPbff3\/q1q2butONiB1Lly5Vsvr999+pc+fO9Oijj9LIkSOVrPAg6wLuvON+Shiy4SPfemAJtvr6EXkmj04KJ32xw5vrgI7iBidHi+nbt6\/CBxeg+OChPuC5kgJKztEtVqxYoW7CnXHGGeq\/\/EF53RMR5dA5dBgEwn2GPn36ZETIk43AJRE714hNPnfKKsKG3SggR0e8zwrN97T1NvKFp+XLl2eQnQO+IfCArIMDKvBAwYKXN9f4rjW+hYcDxXXq1CmNB5Qfl\/bx4DAct49Jg+8y\/nyfHe\/gkg8Cm0Ep5s2bR+PHj6cffvhBfeO9995TBDjyyCOpf\/\/+qm4QC4QPUza6C4VLXRh48PBlKvwZAy4uCLEecD\/RVokVlFgnjFe8UYfUP3xfj97oSgooJgSAaAo41ouLHrgyCAGykupMlsCjvG4JeKS2upRkZc5R3i7sjF9S4DvSrFvdyrMjFs8ppMXB4MEKbhV5Qh+NWJi4aAM\/l28X6hdtrNoAwuAuMwabZcuW0XXXXacGqU2bNtGIESMUSRDZ7\/TTT1dtGjx4sLLy3GcesHItG3xPulByoORILHwJTF5gkgOI\/DMwknJi7KzwluRxIpZeP\/7uSoqnnnqK3njjDQUoJk1oAB6wGhfz8TsUVh8N0QH9FhSTQ79J5UaOIEnBoxWbU1YUtIFHMY5LakcKWDdYREkE9BfRAdk6ytirVrf98F1887777qOxY8e6kkJinUwmlfuE\/+6333504oknEiwcrNoFF1yg5NKsWTNlLeA6gRivvPKKunrLA0EYskEfJclZbzCYyiguciCEjhx88MFpq4o72NKTkKTAAHDNNdcob0SPdetECnnFFQO9vgjjSgqAiIZBCViBuKHwZWEGuVF2EyT+d6fof07XBv2QwmlOwULgUQzKyTGo2JxzaBUWCLtgejmeE4AUbCmsAinbTdR1S+fFUuiyOemkk+iss85SpCgpKVGWAzGn8IB8khSMQS5lg2+wEvNcDRZMzk\/dFNrEUljh7UQK1Al5\/PjjjwqGI444IuNikxEpuGEyrAiU5cEHH1RmHI3iCA5OVgMN1ZWKfXIrcFhwXkmhK5bdHKZ9+\/ZZPi63n02u1ZwC7iSPLtJq6P4qRzBE3+ToZzWnMHGf9DnFIYccQkcddRS9\/fbbaoKNSfcpp5yiSIHQn1BEWC8MRhiBKyoqbC16LmSjDz74O3sLurtalzkFy8IEb55cy1A+0o02cp\/kZIcVRHYAk2+2IiANYqmygusX7PUQi3wx3+3ivQkp5CV\/CT6PCnYrR2y9JEEQCBltlW6dJBaDD58c72PCy9jge3DDrOLmOq0+mZACMZ\/k6hMsAywAXAi0CXO97t27K1JgRG7dujV99tlnSjeLi4upsrJSkQSh+MOQDZOCByn8XS6WWM1FgYPb6pPVXFDH22pCjv01linaIhc45N5TvKNdKzmribb0deM\/1y8EnPYuYlLEpKhf2m7QG7fVz0iQQvp3Qe5CG+BT74vE2HoXcd5JATPG+RPQfC\/pp7x3t2G9EWPrT955J4WMWsd7CJdeemk618JpWHvfsYO20f\/Qv5onaV3jxp572r6qSr1j9S5+81Onl0YgAUzHjh29vBJIWTds8ZGxpaX0SrNmgXwvH5XkAttIkIKPTPDKjQzMi2XHBdOn0z8vvjiNeUXr1vRX69ZU3vZYKvvt\/6jp+vVUTmV05PovM+QyjJ6hOTTcSFbDaQ49Q8MyyvamT+kTOt3o\/Tatd2eUK6Ny+nL9kYS2vrNzpyJ22A9I4YStak8iodq47R\/\/IPp9C5XQXwpP+Syk3jSM5mjYzqUJNNG1SyXr19Olrf\/3v0j20rBdSHM0vO0q60irLN9F3ad16ECfrVnj2g4vBQqCFLJDPOpbdRLEKCraQEVFlUYYeLESKJsiotXUgQ6mtZTQvoBv6w+IgQe\/rVlTszQa5mNCClhitDK79TXtln2w6p8sY9c3J2xM8HB6H1b+119\/NanGuEwkSIHVAJylsnKfjHsSF8xCwMR9imHLRiDvpIgng7lTyxhbf9jmnRRottOWu79uxW8xAjG23nUh0bFjR7jK8RMjEHkEgp472HU4EpYi8tIQDcRqWFjCKSRcct1Wxp2P48hl+6A3KGNSeJRmTAqPgAVUHLh\/9dVX6STz8nBp0Ju\/MSk8Ci0mhUfAAioO3HHn5sILLyTk\/OOc27lYYYskKUaPHp0+es2Yuh0vDwh712oKnRRRxtYJfMZdT1ZpshfjKlStQCRJIVdOcMtPv2rotZNBli90UkQZ25gUBpqKUYBJgQnWlClT6LvvvlMXenAvGZducCmFA13h0ghf8MlVsLb6SIqoYOuXFEFv\/uK0QiSWZFetWkUIlSMfK1IMHTpUkQE3x3DFU5ICN+DwwLLAzE6aNIkmTJhAQUb000mBi+\/5OOxnMKaoIsOGDVO33fQnitj6IUUuNigTEydOjAQpoLxugkMMJAQHgJJbkeLPP\/8kXDnkJxfWwspSgHxRfXr16kW9e\/cuCGz9kALvBL1BWTBzCph4nRRsFQAKIvqxpcnlHKS+uk9RwNaEFGEMPgVLChk7CpfREbVQzikAnls8KT8ANwRS5AvbmBR+NDIC79QXUkQASk9NAO56KBpPFXgoHGlL4aEfoRWNSREa1BkfAu7Yqwp64cSqNzEpPMo4JoVHwAIqDtwR1ILnkQFVa1lNTAqP6Mak8AhYQMWBO5KuzJ07V4UCxYW0XD2RJIUegNlt0vzxxx+r6HgIYcl7F34A43qc9jUKnRRRxpZlZpV1iQ8E5mLvSdeVyJLCNBWTvlTrhwx4x7Se+kCKqGIbk8JBe53yk+k5IbA3gVi2mITx0Q9UjaDDn376qfrKzJkzVZRpxJuVsXE5\/iz7qlwPJxRBfFK5AcjxVwv5PkXUscUyOrI68ckFVpNQLcWcOalI7GgPE9FlrEw854JAxG9EOGefUo7w7D4BSA7tAsXnszHIEoRNPhw\/Rl4HuEn4FvzUm266SeXawI456uEJHR8XueKKK+i1115TpJKkKC8nquWeXyOV0\/dwckZuaEcZW8wVIE++H4HjPJIUoc0pkAktp1IxrHzVKiI++uQ0mslo2SAKCMK7sZIUfJCQd7tx\/l7Wq+e9g3VgUoB4+nERkOaXX35R\/66TIg9xzgxRxdknRPr+u3iUseUw+XZzitBWnz75JBqk0EczN7+XT3biSMe7776bHuEx0cbjRAocHJSHCZ0sBasTWxTdUuD3qFsKec7SJHVuvrDlBDM4DW3lPjXofQorEw\/fnpPE8FyALQWOiyP1GEBDEhM3UmDugYSJeHBgrl27dirHw4wZM1Q9nFiGcx7IZDWYdxT6nAIpyJBbg58oYYu2IPvunXfemZHzO97RNnYOwi9Y6KtP4SMWzBfDxD2SS7LBwJibWsIUTm56UJi1hol7TAqPOhKmcDw2rV4XDxP3wElhlyGVJaanquV\/N0nWIlcl7LKN1kUz5G00u\/rrIhw9zxy31SkzLMroueHsyhc6tk6yc8PdL7ZW3wyUFLJhdsd8\/SqzzGkt16\/rQgL5rp5R1a5eN+E4tcdvH6ySUGJDEldx5VPo2NaVFDLfdl30whMp2ApgVB84cKAKEoDVGkQM59S53Bi\/pJDE6ty5M82bN4+w6carRagfdWMTjjOevv766\/TAAw\/QqFGjVI5mPOPHj6dvvvlGhcqRl430lS30pU+fPukgW3gXIzEemVFVf4\/7x\/8+ZMgQQgIRrIzZjeROpEB\/sCoEwh944IEq8BeWjoHvyJEj1f4Ib2AOHjyY1q5dS4lEIqNv5513nsKqZ8+e9MEHH2TdPyhkbJHRFf31g62eyxsyw77T2LFjVX2MLfZJgJEnUkBROH0v7v1CMVg5+M40J6L3QwpWsIsvvlgpJP4OBcD6NXY5cWDvxRdfVIqj58bGUqkkK+fDbtWqlVI2LNdyrm8OuShTEaNvKMfpjq1GZ\/zG3+FRid\/r2rWrIqRMp6xbNCdSQHAYaBDAoX\/\/\/vTOO++oZUnkKZeKAByws37jjTeqvRkMTHjwbQwcGAgg5L59+2YMloWALfacJEYSW8gT8rWzBm7YAqfly5enUxbzAC+x5To8k0KOmFYByng0ciKF3C1Gx\/l8EYOANfQzzzxTXTHFA7CYFM8++6zaeYby4UFdSFAPSwGrgo28xx9\/XP2GDSBkRYLVQP5oKBFOWV5++eVKqTgPNke7sCMFrCALo1+\/funwjXBfUC+TDiON05zKzu\/lnN3XXnstLV26VCWBB7ZcN2PJ1qSkpEQpz7HHHpvxPT7fJdlQSNhiEJGpfCW2kDNwAPGtdMsNWzl44IiQLieZMjqR6t07+5jHJ59ku2Sn\/53mavrWrfTE9u00v7SUemjn2hds3063bd2a3fDa9\/HukooKerq0lFo2avT3d3AWoaxMWQBJmm7dutHdd9+tFB2WAgqIDuD54osv1MYbDpDhaMZjjz1GixYtoqqqKqW4cOu2bt1KF110EQ0aNIiQjB0jbzKZVCPGihUr0sruZClQN7tS2PjjmKZpUpx9Nl3fogVd0qIFLdm9my7bsMEam5Ejs0e6Wlx2J5N066ZN9H1lJa2srlbv4+G6WhUV0VAiWo+UVqedRs8880zaxVOKUlpK0\/\/zH\/qpspLebNOGDm3SpAZbHBWojZQSdWxxFs2KFGyh3UhhZ0WscqTbJbCHTBOpsrJsUuAgkv7UHvBZkkrRZUVF6tfjUyl6OpWilom\/k10tSKXotqKibFLUvj8dvjpR1nsEItaeR4A1gjV46aWXqLS0VJ1+hT8JEoAUIApGzcWLFytXCmWee+45ZSm+\/\/77tEsBSwGXCS7YMcccQ3v27FF+v3T5dDfIyn1ytRR9+9L1qRRdkkgQ4zM\/maQeAhfgtWDatGxS1OICmo9IJunL2oEC2F6TTNLIRo3o8epqml1URN8mEmq3F7v2cJcyRrxkkqYnEvRTKkWwoYfyt0EKcfgpytjaWYowSJFpKVLmZ5+YxYjQh5F3xIgRWRMfE\/dJTmAl95i9OG6BERFxnPBvcCUmT56sTrY6WQqQAi4FRkQmDCbAcJ9Q13HHHUdff\/21IgpPrjdu3JgOFlaXOQWTyc19shvN4Mqdc845tPfee6u5BFwp9Bu+tPwvBghYUli69u3bZ80pMNGGawnyFBq2dnOKMEihBq0FC9Sg5WlOoSu8nKjypNKEFPqcAg3iCY\/uG0p\/G5PIbdu2KeXm809yTgH3CasUKMdzChAAyr5r1y5FYkT04+9jlebzzz9XlgMrVBxyE0oI68MuGB8x5zND8o6FvpTrRgpM1PUH30OfcBAOLh6sJE7+vvzyy7Rp0yY67LDDaOXKlRmvYSyD9ZTzOqw+MS6ycH3AFri5uU922PL9GGCSZWHnz89Y2va1+pTtV0X\/X\/jUp37y0k\/L67JP4ed7UX8nSGyd+hom7p4sRdQFZNU+OXHj8\/p16UeYwqlLO\/ndoLP8yDYFjW1MCoFALgUXhGLJOgqJFLkIPhw0nqb1hYl73i1FoQkuTOGYKoxduVxk+alrm\/y+Hybuec+OKhe\/MHHkv+PP8imtrqYNtcuV7auq1E\/dKyrolWbN1J9bVVfT5trfW1dVUd9du+jZFi2MZID61jVunC7bPJmkHbXLznoFKItvfdu0KTVKpahNdXXGu1blq9q3V6tpYT86lvg7\/l8k+tY02YT2q2xJp1T9ktWP1VRGB1M5VVW1o0WND8tofhmtVu\/gwW\/8PvDB3\/Ee47qQemV1He\/jQTm7h+sqp5rblPzgXX5v6qJFgadCyLulcEvPdNpBB9Fna9ca69On1Pu\/\/6uJ4uH3KacyKrMQ1lLqTt1pqVG1so4lJSXUo3bD0ejlgAq5YYvPrCsupub\/PJW23D2JEge2oqL1G2hPq\/2paOMm2tOmHZVs30SpqiT9eUD7dKswXjXasoVK1q+lRPO9qaq4hCiZokbJPbSnuogaFycoUZSg4vW\/0S03jaNx\/\/pCvZtI1GyJ4f2SDf+moj82UaJlc0o1bUpFyWqqbtyEaHcFNfrjD2q6bTPtOLYrFf\/1J+1oU7O8zO8XpZJUvGs7TRs4gD7eu2ngNyEjTwqYTTmSwzosLSmxVJueu3dTZa2F4TIX7NypRkD8\/cSKCsIeOv\/G9R5dWUnfNWmivrOHiMpqLRF\/RH5Pfn\/gzp30R+PG9O\/GjaldVRX9XFxMh+\/Zk1E\/\/h3v5+MKq4n7dGnbtrZ4BsTNnFcTNLaRIEXQ6ZlyLoUC+UChzdeiAmveSRELLreqEHSWn9y2Nhq1550UgCEWXDSUIW5FDQJ5X32KBREjYIpA0HMHu+9GwlKYgpKvcmGukeerj1H\/LsuAT7PyRTHd0zC56+\/W15gUbggRUUwKA5ByXIQDLCOvIQfKxt2HXMxJY1IYCDMmhQFIOS4CGeC0L4JjjxkzhnCODaQwWXb22rTIkAJHsMeNG6dCV8pj6OhQUCmATXNQ6CAWMimijKsXZWUZ8CFESQqOMI\/6cPwfV5DrojORIgWuinbo0CGdvgln24MkhRchyLKFToqo4upFHg2WFIg0jrAuTAQmBYfT4YDH+sV1mNDZs2erC0N4eLLlFLZfhs3hEDjy32SUh0InRVRxDYoUQW\/+JvKZn0JehOUQ8ZhIcSh2BBbgBzfmYDI5iYpMHSvP+CCyB0Lrw3zOmjVL3cTDzTkAxz4oonqAWIjMwa6aXq+0UjopEDNg+HAvIg237MSJ6VgF6ZwcUcTVCyp2liInE+1hw8zvaHvphElZu2QinAsC94xxZxkP55Owuuklw10CJETfOPvss1VwAwRIwwMlr6ysVGFtQIrNmzera6r4OywF7kWPHj06I0Q9WwudFMhHMXeuSQ\/zU6ZXr5pkLXhkPoqo4eoFHTtSoI6gN38jNaeQiVrY9eHsNW6WgpO0SFI4WQqZARXfQmADLPVZJS+vD+6THBw4oBsUKp+4+iGFl3f8lo0sKXiTBisJJnMKnRQIWICoFnyZXU8FJlN4Wc0p5NykPpEiSrh6UdowZRAZUngBKOyyYQok7L4VyvcgA7uok0H3ISaFAaIxKQxAynERyKBB57zLMb6eq49J4RmywF+ADELLjpryECEw8J4WSIUxKfIvKMggtDzaUSGFzFJkFRAXK0RYKSkrK1N7EAi23KVLF5UXg5dr6yI6pyMghUyKKOPqRV58IBBR461WCL3U5VY2MnMKuVmG9XQIEw+WEjkfBdL7YgNOKnBQpHACqpBJEWVc3ZRT\/t4gSSE3mXgzDqDACnTq1EnllRgwYIAKrY8HS62YeCHgMDbgOD6svkIBpUDOCj4CIiOO2y3X6kdAsH8R1gUXL4piUjbKuCJVgt3RGr1v4ZJi4sSauCO1OQwyGjNpkgnuNWX8vC\/ekaM\/9hBgFfBgVxrJO+BPggDYj8C+BQIQY2caloI3oEAmlJfhMdnt4oQqXN\/UqVONj4BYWgpTbPzgIlG3eh+\/O30f2CEEP1GGVY0ark5Ha6xIEd6cgqiGFFanPbSAZI4M8fO+9g5cJig+jgLDb6yR\/STlPuGgIKyGFSl4TiFHRW4rpx2D28W\/48yTvnvudATE0lKYYuMHFwm03Skcp+\/Lw09EyhWNIq5OR2usSBHe6tOqVTWkqE2YktGYcvvobVkE8fO+9g6fdkVSEh7tIdBly5aps0k4rxQUKZwshX4EBG5Wlvtkio0fXCS4Vu\/jd6fvFwiuTkdrrEjRIPcp+AIJEq3wJRGM9EjAgrwCGLGZFJyZFUCBRNIS8DkfAGtlKfC707Fy\/QhIIc8pgEFUcdXnFGir3R3reEfbfCYTSslCXn0KBaAQPhKmDCKzJBsCrr4\/EaZAfDeynr8YpgxiUhgoU5gCMWhOgywSpgziYGgNUsUKs9Nh7RX9P7jf0jDpiuJMAAAAAElFTkSuQmCC","height":0,"width":0}}
%---
%[output:352ed959]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAMUAAAB2CAYAAAB4ZlcrAAAAAXNSR0IArs4c6QAAHnNJREFUeF7tXQlwVNW23R0TZApGhkgEMSiDYJUTqIAyaTGJgrNI1ZP4RUVBcQJnQ1RQEXAupBQD9ZzLocQB1K8ofhEK8amlXxyJCAQR8YnwQQzp7zph99t9cu7UffvmdueeX798pM89w9pnnb3PtHcsHo\/HKUoRAhECCQRiESmi0RAhkIxARIpoREQIaAhEpIiGRIRARIpoDEQI2CMQaYpohEQI+Kkp7rnnHpo3b16iyCOPPJLmz59PrVu3TvwNeZCuv\/76nAX\/ueeeo2effbZe37nDW7dupYsuuoj69OnTIDjocnr66adVW5wS+nXjjTfSXXfdReedd55Tdk+\/f\/vttzR16lSaOXMmde3a1dO3MjPauHz5crr77rupWbNmKZcjP0xJU7CQO3XqlNSYFStW0NixY0mC3phIcfjhh9PNN99cTzgPP\/wwvfPOOw1CCuC\/cePGhJwwGC+88EKaNWuWLTFYxq1ataJt27ZZEj7VUegXKVKt3+67lEhhN9BBjOuuu44qKyvVDNBYSIHZ6sADD6QzzzwzaebD4Hrsscfojz\/+oMLCwkA1BQ9saGmpGdzIBHJEPnwLeTqRyOvgzClSOHVm586ddMMNN1C\/fv2UygWwGBD4\/1dffVVhd+mllyYGB+fn3zCwmFDIq\/8uTTS0ZcaMGdSrVy+aM2cOjRw5UpV\/4oknJql7FjCbdtKc0OvDby1atFAz+2effaa03vvvv29rArIKHz16NP3yyy9JdeM3pKqqqqQydJNGtsMJE9bIPBAlnm4GpxtScJ4rr7xSyRPtkyYw\/4762IQ+7bTTkiwHuz7KcfTSSy\/Vw1eaRcAD5ifkgSTHgG4+scnHOKRi+nnWFPoAMwlBNvTBBx9UoLFJxbPXmDFjEqSRaw4TGJwX+aT9\/uuvvypTAMJggZnsezsBmzQbCCqJ6TTQuM1sI5eXl6t1FYQ5ffp0GjduHEnB64JkEqAe2MbAzAqT9evXJ9niVtrAqs06\/qZ8uollkjkPeCu5OvVR9gNyRHk8acmJddSoUUmTLNorSS3rAWlkOU4TuBVGKZHimWeesV3Y6KSQAsb\/ZpBha2Ohpc9C3Fjk0+uSgB1zzDH17GN9kMh\/t2nTRmmW2bNnW24G6Da4EyGYqLzYw4AeOHCgMlcglIULF6p1hj7Q9XIZs9tuu41uv\/12W0ykeeqmfZxHJ5\/VwlQf0CYimXBy0kCyXEkKyAWagM08OZj5Nzkxyj7LMhctWmS74eEWq5RIIdnoRlOUlpYmmRSy05glsDhH0k0ZXRXKuqAWQQrTDoYUjpzlvvnmm0RdervZBHESrFN\/MVvB3IKA2XRiM9I0OXDf8RubHyjDChPdtHJrOvF369ats1006+Yv91cngQkn0990U4\/7qGs8q9kfxJXjQN\/htDOzUjGd0F\/PpHBSSaY1hR0peDtOCpvJAZND7pzoA9KqLZIIWORy\/U5bp7pqdjuz6IKpqKhQZJ07d64ynfQNBxaynARMW4smTOT2pRxwduRwSwjW4pKoEgPZXidSOPVRJwWvD++8805lPfCaVNbPGguTBpPj7bffNm7JSiJ5JYdnUjgNHDe7T3brEkkq1GW3\/29FCgYP9ihUKu+F620zDfp0NQVmNpRRXFxM33\/\/fWKLVl\/X6EK3M9usZm9pZlppb7cmk5VG4L\/rZy12pODFuV0fdVJwO4866qgkmZlkJE3iH3\/80facIpVzjJRI4fWcwmqhzYsouaYw2ZPy0EsuAmFvWh0A8Uwhd0RMA8Rp4e9GW1jtgMgZym6xL9tqWlPo5qa+prAjspc1ktOiHWVhYsGCGBpYNwfd9hGbCTopeG2Gw0IpM9N6Rk6qUlPoawqnycRKtimRggvTbX6rE219S1YOFqkSuVx5+Gf3u50px+S54oorktYzuk2OOk2D17T9aHUqr5PC1C45cPU+YRBgOxeEwK4XLy55CxJtlJjouOtboYwjYwATVE+mb5zMSzkhmbapvfQR7dEnNCuZ6f2wMzv1bWC3ay6JT1qkcDOLRnkiBNwi4LRedVtOuvkiUqSLYPS9bwhAU+GQs6HvyUWk8E2kUUGpIsDmJL7XL5SmWmY630WkSAe96NucRCAiRU6KNepUOghEpEgHvejbnEQgIkVOijXqVDoIRKRIB73o25xEICJFToo16lQ6CESkSAe96NucRCAiRU6KNepUOghEpEgHvejbnEQgiRTyopnukQM3XfVHQJlARL\/8ZboLb7okiLboF\/v4Rqd0ueOlzfJWaKpleKnPLm8qssHLxmHDhqXlQka2yUk2Vu82rG5Vy7JNnmC8YqeXgfa++eabNGnSJE9FWZKCbxc6PXTxVJtDZv1qN26Jmp5emvwoMSBeH5TYNSmspHAjG5DooYce8vTW3A4Lt7Jh8spJ1Y1s\/CCFbH86vraMpOjfvz8VFRWpd9i49w7nAIceeqh6NGO62iyv51o9HdRndzeD1+oqsVWH9fv+UlPIK8V8bRogwlMF\/Boh4Tq0\/rrMqgx5Td7uiajEg+tl7xTwE8VeTlDexRdfnJjVJD4oA5fk8vLyqG3btnTcccepd9zXXnstrVy5UsmqurqaunfvTo888ghNmDBByQoJkRbwzh1vUoKQDctMdybBBNWvxTN5dFLYjRcrvLkMjFG82mQPMUOHDlX44NETXza0m\/CMpMAgZ48Wn3zyiXr9dvLJJ6v\/coXyiScehqBz6DAIhDcMQ4YMSfKKJxuBhyFuZjEeyLpHOytSyNkRoPCA5rfZehv5kdPq1auTyM5O3uBsQJbBThR4omDBy4c8\/L4adSGxc7hu3bol8MDgx0N9JFyA4\/YxaVAv489v2PENHva0a9dOeQl56qmn6NZbb6WvvvpK1bF48WJFgB49etCIESNU2SAWCB+kbHQTCi8RMfEg8QMq\/G9MuHgUxOOA+4m2SqwwiHXCeMUbZcjxh\/rtPDYaSQFbFAKABwVc5cUjFTwThADZ9YvOZAk88uuagAe46SGSSW3zrGbKnyopUI8kmOklnhWxeE0hNQ4mDx7gJm8T+mzEwsTjGti5\/KJQ14imNoAweL+MyWbVqlU0ceJENUlt2bKFxo8fr0gCb36DBw9WbTr99NOVluc+84SVadmgPmlCyYmSva\/wwy\/5aElOIPJ\/AyMpJ8bOhLckjx2x9PL18WckxeOPP06vvPKKygvvD2gAEliNx\/j4HTOZPhuiA\/rLJyaHvkizI4d8kG\/ye+qVFDxbsTrlgYI28CzGvkitSAHtBo0oiYD+wiMga0fpb9X0wg\/1os57772XpkyZ4kgKiXVtba0yn\/Df\/fffn3r37k3QcNBqZ511lpJLy5YtlbaA6QRivPjiiwSsGMMgZIM+SpLzuMFkKj23yIGIMXLwwQcntCreXUtLQpICE8Bll12mrBHdv60dKeSzVkz0dpswRlIARDQMg4AHEDcUtizUIDfKaoHEf7fz+Gd6KijJY7XuSGVNwUKQGoj9TrE6Z3cqLBA2wfR8vCYAKVhTmJwnW9mtevu9aApdNn379qXhw4crUjRt2lRpDviZQgL5JCkYg0zKBnXwIOa1GjSYXJ86DWg3msKEtx0pUCbksWbNGgXDYYcdZvmYyZIU3DDpSgSD5YEHHlBqHI1irw12WgMN1QcV2+R276ft3taaSKEPLKs1TMeOHevZuNx+VrmmNQXMSZ5dpNbQ7VX2Woi+ydnPtKZwYz7pa4pDDjmEevbsSa+99ppaYGPRfcIJJyhSwN0nBiK0FyYjTCq7du2y1OiZkI0++eDfPLnp5mo6awqWhRu8eXHtZIFw2y1JoQ8Q2QEsvlmLgDTwn8oDXH9Ur7tV5Mf4pofzVg\/tdY2RyjmFNB1Ye0mCwPkx2irNOkksBh82Ob7HgpfJAzBhhpl85drtPrkhBfw8yd0naAZoAJgQaBPWescff7wiBWbk9u3b07Jly5R8CwoKaPfu3Yok8IoYhGx4YLEs8W\/pgtS0FgUOTrtPprWgjrdpQY7zNZYp2iI3OKzOnhr1ibZXn0j6giz6d3Yh4PbsIiLF3u1CP4N+ZNdQaRyt9bL7GSgppE3n5vCucYgr\/V5GuKaPoSwhMFJAdeEE9qabblL1m7x\/+9u1xlFahKv\/cg6MFNLVIZ8bnH\/++YkIO1OKi+nFli3972EISkQAmM6dO2ekJU64DujUiY7f+Se9WNgiI\/U3dKGZwDZQUvA1Cd6tkQ541xcUUMeaGvq9zwBqEt9NMYob8S6jBVRFpUm\/PblhCJVSXaQgp9SZ1qosmzZtUjs1SJVURoOoLlqRnmo6dEj86X82dKFbOjyelGUQvUfT\/v4\/uzTv669p8m+\/OTUtpd9BCjtcsUO27qefaCY1peNoV1IdLdt3pRqKUY\/2zWkB\/YMW0H\/BEX0iD\/qF\/rlJ59MzVE0ltGbNV9SjR11gx3Hxf9I\/6EklyZpYPjWJ\/0W7StpT8+qNtKukhOKbtlIN5dNWKqIyqtxbN+qP0YG0kZ6MX0BxilEe1Rqb8NWaNTS1XTtatm6dmya6zhMaUkB4tHewmwY4iIC\/64RAT90SAnn9\/t5N3ahz3bq6bVK\/kxMpcK6ByYYT46i3w29cgpAL2pyfv55++OEHX2ENlBTsLt5kPvnaq0ZUmJP51Iig8K2rgZEiWhD6JrOkgiJc\/cc1MFKg6W6P2f3vZm6XGOHqr3xjnTt3Nq9o\/a0nKi1CIG0E\/F47WDUoUE2RNiohKAAL16CEE4LuhqYJjDtfzZHb+X4fXkak8Cj2iBQeAfMpO3D\/+OOPE0Hm5UVTvw+FI1J4FFpECo+A+ZQduOM5wdlnn62cWXDM7UzsvoWSFJMnT05cw2ZMreK6+YS562KynRRhxtZOCIy7HqzS6ZzGtWBFxlCSgtuHDuPFn\/7sMJWO+vVNtpMizNhGpHAxSiUpsMCaPn06ffnll+pxD94o4wEOHqiw0ys8IOHHPply3JaLpAgLtqmSwu9DYVw0CcWW7Nq1awluc2QykWLcuHGKDHhFhueekhR4DYcEzQI1W1FRQeXl5eSndz+dFHgEn6nLfi7mDccsZWVl6uWbnsKIbSqkyMThZWzatGmhIAUGr5Pg4A8JjgIwyE2k2LFjB+H5IadMaAuTpgD5wpoGDhxIgwYNygpsUyEFvvH78DJr1hRQ8TopWCsAFHj3Y02TyTVIrppPYcDWDSmCmHyylhTSyQEepsODoVxTADy3jte8AN0YSNFQ2Eak8DISQ5Q3V0gRIkhdNSVI3EOtKVyhFXCmIIUTcNdCXV2QuEek8DgUghSOx6bldPYgcY9I4XEoBSkcj03L6exB4h5KUpg8Bdotmt99913lKQ\/uLPnsIpURwuXYnWsEKZxU+uD0TZix5babIjAFiXtoSeE2LJO+Ves0KKx+d1tOkMJJtS9233kJeeUWE6d2ei0nIoUBUTvB6fEhcDYBv7a4MMhXP1AkHBC\/916dJ4pZs2Ypj9PwYyuDgrBfWzhmQ+JyOLgIfJXKA0D2xZrN7ynCji0sAkR44psLPDyCnIxilZXxUJxol5X9hx0mFc9xIeD9G97O4fwASc5CbD7h7+z2BQOf78YgYhAO+XD9GDEeYCahroULF9LVV1+t4m7gxBzlIMnrIhdccAG9\/PLLilSSFFVVRHu55zRhNsjvuDkjD7TDjC1CCECe\/D4C13kahBSIitYg0tIqXbuWiK8+2c1m0nM2iAKC8GmsJAXfruXTbty\/l+XqMfCgHZgUIJ5+XQSk+e6779TfdVJkyM+ZL2LBZCOvPoUZW3aZ3+Dm09Kl4SCFPps5rSn4Zidm8zfeeCMxw2OhjWRHClwclJcJ7TQFj0zWKLqmwO9h1xTynqWbNUVDYcvBZnAbukHNpzjiQYUsmVQ8bHsOGMNrAdYUuC6OMGRYVyCgiRMpsPZA8EQkXJjr0KGDivcwc+ZMVQ4HmeH4BzJwDdYd2b6mQDgyxNngFCZs0RZE4r3llluS4n8HuqYIIylCxtGk5gQpnDDjEHTbgsQ9lFuyQQPupb4gheOlXbmeN0jcI1J4HE1BCsdj03I6e5C4p0UKPX6ZVSAWPVQtS89N4Ba5E2EVbTSd0SBfoLkp36tw9Dhz3Fa7QJeyP4wxx9fjrWjOk0vY2snRhHu62FrVlxYprCKQyv1lVOxmsJkaKGNa62WmQwS5o2SKge1VOHb50+0DC97Km0kuYesV93SxTZkU\/NQPs\/qoUaOUYwDs0MyfPz\/p7bNdA50EJxmPnR5szaEeucuDQ7Wvv\/6afv\/9d1U3ziQQwhgPi+644w7VP3wDxwbYNZJ3pfTdLPRlyJAhCcda+BYzNxKHouWDPblTg+1Y7Fphuxh\/HzNmDCFoCHbDrGZ+O1xwkIhyQPgDDjhAOfvCdjHwnTBhAj366KOEqK3bt2+nI444gj7\/\/HPVRtm3M844g3Ao2b9\/f3rrrbcSgeRZ4NmELb+3N2Grnw+hf3bY6rG8MQZw1jRlyhQlK5THMtMxcqUpMKgx0PDWFwNDhm91ownsSCHjX\/NA7d69O1VXV9Off\/5JP\/30Ey1YsIAQpnjx4sVJpMD2qCQrx8Nu06aNGmyYXTnWN7tZ5L7wY36pKUwhgmUoZJBww4YNavDhu169eqmg6TKcsq7RnASHiQZOG0aMGEGvv\/662opEnHKQffPmzYoUX3zxhRr4KKtLly6JOOCoGyT99NNPlZCHDh2aNPllG7ZSJjq2wB04yeSELSbJ1atXJ0IW8wQvJzBZBoeQdkUKOdNaqXE9hrFsvMnu5TtFyMezMa8x+IAJIODm6vPPP58gxW+\/\/UYlJSW0ZMkSisViKqb0SSedpAZRbW0tHX300eraBhJuzuKUFDMyYknjQIgjGHEeK1KwJuIY0GzbYwAyKYAFyrfru5XdyzG7L7\/8clq5cqUKAo\/y0E60CeR44oknaM6cOXTNNdeo\/oB80GCyPr7TJfHOVmyhKbj\/OrY4TjORguO5y\/4ztiAFEke+1eVkChmNPLH4oEH1D++WLq1nbt1z1FE0b9s2erq4mPrsvXckM23ds4cu2ryZ+jRtStd\/+mnS92oGrqyk+cXF1HqffeqX\/fdAR9mcrrrqKsrLy6P99ttPrUfuv\/9+FUgdmgImFUjQo0cPdZgHrQJzC2oRgdSPPfZYuuSSS2jixInK\/Q3uOZ1zzjnUtGlTNbg++ugjNQidNAWuekhTykQK1iIr+vShsZs3G7F5bts2emjbNqps1466Ll+e1HcllN696X9376bv9+xR3yOd\/\/PP1KWggMYXFlL\/Zs1oZHU1tc3Lo0Wff67ufOmk+Oazz2hRSQkd2qRJfWxbtqR5X36ZFdhKUiSwXbFCaUMrUkB7y+D13FGrAY+ypB9aTKaol6+Y1JGitLQ+KTQ1JV2IHBmP0\/x4nFrH\/hMbDQ35Nh6nC2MxuiIep\/NwS04kRYpHHzV+l8i2dCmt2LRJAVBcXEzDhw9XA\/++++5Tt1+hTjHjgxSYRaExcOcJmgEzK0wOaIrRo0cn1hisKVAm\/I\/iFBdlQ\/M4kcKNpkgIrrSUxubl0dO1tdRHw+W5eJweisUI9y67arhAcOMPO4w+2jtRANvLamvp0rw8pQX1dMoppyjNiDUMCxfYfvOvf9EiIjrU8I26DVhZmSBSmLENmhTWmsLhmofcEsSAGz9+vFqgsK3O6snOrrZbU+gqDXlx5Xvw4MHqfhJmgnPPPVdpiwEDBthqCpBCDvaavwNLYnEFErVq1UotVmGaQes4kQJ5pGnFtqdcU+izmb7WAjZ2di9Ieuqpp1Lz5s1VG2FKwWzgtREGCQSH9ZDdmgK\/AX9cj5Ap27A1yYT7kAlNocvH9ZqCbWIWuGmhyvdoTIMCFVvtpZtW\/8iPO0iw\/YcNG5ZYb2BAFxUVKQ1htaaA+QSzCwtTJGgKXPfGBTPsWCFhl+aDDz5QKhS7OOxmE4MRsyh2gHhnTa6l2E5Hm3j3yS0pTHYv6kOf0La+ffuqG7rQfC+88AJt2bJFXX1nUqCNuODIu09yXYfdJ74LJgmRjdhik0JuDuCypxMprLDlNzF2awrGS677gK2rhXY9PZ5lf+Bbn\/rNy1S64fXwLpU6sukbP7G163eQuAdKCr8jzrgZPGz+ycWUm++s8gQpHLftbAhc0Ta\/sW10pMiEI1y3g8bPfGEjRa7g6iSjIHEPTFNkIuKME5CZ+D1I4bhpf67g6tTXIHEPLDoqb3LxViP+jf\/HwpjTWdu30\/r8fOpYU6P+y0n++30aWA+\/gfS+EVO9nB+plKqo7mWeTPJ7fIMk68e\/8fcP87vU+76UflT5D6bkbWhZ94b8\/HrlOQ0Ct7874VpT05HO27XGsX4TNuib3i+rdpnkYvV9nGIU0yJAWMkV7UIb+L96\/TM+\/ND3UAiBagp2JoCOYUelX79+iShFK5o1o95t21LVffNo3wMK6a+i1pS3cRPVtm1DTTZtpJrCItpd1Ib2FOxbTy6F67+leEEB7W7Vmgp2bqfY3oG9u1kh1eYXqN\/22fl\/lBeLqzxIw4cPoyVL3lQv1Fus\/552tSmheH4+xfbsobzaGqJ98ii2Ywc12fE7USxGfxW1VXXX7FvnMIFT8+\/WUDwvj2o6dKC8XTspf\/s2qv1zD\/1VUkK1BftS01+r6ZKxZbRkww9ux7mnfE7hrdQMu3YtvXb\/w7SrSQvquuy\/qfDnair8eRP9m4qoiP6t6nuPBtHEX2Ym1V1GC6mcprlqTwVNowU0Tu0OwikEEqagpTTY1feDaSlVUXJ8EnyLMqwS6hpbUkLL1q1zVYfbTIGSwi7izIBOnYwztNuOhD1fpp6wOplPIEX7Xe1oU9Nfwg5Ryu3zG9vASNFYFoQpSzbFDyNcUwTO5rPASIE2+B1xxn84srPECFd\/5RYoKfxtelRahEBmEAhs9ykzzY9KbUwI+L12sMIu0hQuRlWQe+QumtMos7AM+GYrPxrTzXI37\/6dAIxI4YQQEUWkcAFShrNABrisifcP7CgbV3cysdEQkcKFMCNSuAApw1kgA9z8xaMxvI3BoyCQwmlLOpVmhYYUuDI8depUdW2c3znjSi+SXyGAvcZJYECzmRRhxtXLgGUZ8CVESQq7Q2EvdXDeUJECr+cOOugg9bwUzy79JkUqAOGbbCdFWHH1Io9GSwo83oGLF9YOTAp2rcMOj\/XHTFChc+fOVY+HkHixZee2H6\/V2Mkyu4yRf5MeH7KdFGHF1S9S2N2U8FJHQlM0ZHwK6e+cPXhgIcWu2OHWhlNVVZWyI6E+KyoqqLy8POF3St7\/Wb9+vfLmAZMLnjBmz56t8gE4tkHhQwrEwis6NtX0cqWW0kmxYAE8kKQCdzDfTJtGVF5eV1eYcfWChpWmyMhCu6ys4VzxWwUT4VgQeHOM98tIHE\/C9NJLur4ESPDEgUf+cHQwadIk9T0GOZ6r4uksSMFvsPFvaAq8kZ48eXKSi3rWFjopEI9irxcdL3INLO\/AgUQcGUrGowgbrl4AsSIFyvD7RD9UawoZqIVNH45H56QpOEiLJIWdppARUFHXjh071Faf1EC5stAOK66pkMLLN6nmDS0p+JAG18vdrCl0UsDbCDxc8MN2PRSYDOFlWlPItUkurClYY4YJVy+DNkgZhIYUXgAKOm+QAgm6b9lSX5AyiEjhYlQEKRAXzWmUWYKUQUQKF0MsSIG4aE6jzBKkDCJSuBhiQQrERXMaZZYgZRAaUsiIRSYfn9ghwg5UaWmpOoOAN3J4AOTY2TiDSCfZXQEJUiDp9MH0bZhx9dLXIGUQGlLIwzLsp0OYSNg1wQEbDt8Q3pddSeoB5dMlhZ2AghSIl4HiJm+YcXXT\/obYFg8NKeQhEx\/GARAc2nXr1k2dYo8cOVJ5i0DCViv8fsIjOQ7gELQESb8CgkGB4PN8BYR\/93IFBOcXQT1w8TJQ3OQNM67w2Wt1tUbvW5ATUyw+bVqdK36+FyBbU1HhBve6PKl8L76R5gvOEKAVkDjgCg6gQACcR+DcQmoKPtgDmZCfYw3geza7OLgKlzdjxgzXV0CMAnGLTSq4SNRN3+N3u\/qBHVzwEymv5YxV2HC1u1rTsKSgvV6pTLc9TPEOrGiSyvfaNzCZMPBxFRgny3Wyr1DmEy4KQmuYSMFXQOSsyM0EKfTf2XO42ysgRk3hFptUcJEYW93CsatfXn4iUqZoGHG1u1rTsKRYu7ZOU5QmO6JSf9OCjNiqjVS+177h2649e\/ZMzPYQ6KpVq9TdJNxX8osUdppCvwICU62e+eQWm1RwkUCbvneSTZbgane1pmFJ4RC0xb39lH5OfkCCqKP8sAgz\/fLly1XcMszYTAqO0op1BUjEkTXlPR82n0yaxG5NoV8ByeY1BTAIK676mgJttXpjHeyaIkSkSJ9WmSkhSIFkpgfZX2qQMgjN7lOYxRakQMKMQ0O2LUgZRKRwIekgBeKiOY0yS5AyiJyhNcohlp2dDuqs6P8B33\/kf45Dsd4AAAAASUVORK5CYII=","height":0,"width":0}}
%---
