# TP1 : Séries de Fourier


## Première manipulation de séries trigonométriques

__1.__ Afficher les fonctions sinus et cosinus entre -10 et 10. Combien de périodes observez-vous ?
__2.__ Quelle est le décalage temporel entre les deux fonctions ? Comment reconnaître cos de sin par lecture graphique ?
__3.__ Afficher à nouveau les fonctions sinus et cosinus en changeant la valeur de la fréquence à 440 Hz. Quelle est la période ? Zoomer sur une période, que devient le décalage entre cos et sin ?


On fixera à nouveau f et T à 1. Les coefficients a_n et b_n sont tous fixé à 1.

__4.__ Afficher la série de Fourier correspondante pour N=1,2,...,100, 10^5. La période du signal change-t-elle avec le nombre de terme dans la série ? Selon vous, est-ce que cette série converge ? 

## Analyse et synthèse de signaux péridodiques par Série de Fourier

__5.__ Afficher la fonction carré entre -10 et 10 et une période de 1. 
__6.__ Calculer à l'aide de la fonction dédiée les 100 premiers coefficients de sa série de Fourier. Rappeler la formule des coefficients a_n et b_n pour un signal carré et vérifier que vos résultats sont en accord. 
__7.__ Calculer les coefficients du spectre du signal carré et l'afficher. 
__8.__ Afficher le signal carré reconstruit à partir de sa série de Fourier pour N=2,5,10,100 coefficients. Que remarkez-vous ?
__9.__ Reconstruire à nouveau le signal mais en utilisant seulement les coefficients de rang 50 à 100. Que remarkez-vous ?

## ...

__10.__ Afficher le signal *signal_exotique* entre -10 et 10. Est-ce un signal périodique ? Respecte-t-il les conditions de Dirichlet ?  <!-- to CLAUDE : Yes it will, a Weierstrass like function would be very nice no ?>
__11.__ Appliquer la fonction de calcul des coefficients de Série de Fourier et reconstruire le signal pour plusieurs valeurs de N. Que remarkez-vous ?